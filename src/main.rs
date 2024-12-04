use std::{
    collections::{HashMap, HashSet},
    thread,
};

use asts::SubreadsAndSmc;
use clap::{self, Parser};
use gskits::{
    gsbam::{
        bam_record_ext::{BamReader, BamRecordExt, BamWriter},
        get_last_pg_from_bam_header,
        plp_counts_from_records::get_base_idx,
    },
    pbar::{self, DEFAULT_INTERVAL},
    phreq::quality_2_phreq,
    utils::command_line_str,
};
use ndarray::{Array2, Axis};
use rust_htslib::bam::{
    header::HeaderRecord,
    record::{Aux, AuxArray},
    Header, Read,
};
use time;
use tracing;
use tracing_subscriber;
use utils::build_target2idx_and_idx2target;

mod utils;

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(short = 'q', help = "subreads.bam")]
    pub sbr: String,

    #[arg(short = 't', help = "smc.bam please set tn-delim and ch-idx")]
    pub smc: String,

    #[arg(
        long = "rep-thr",
        default_value_t = 3,
        help = "only insert a base when the repeats cnt exceeds this value"
    )]
    pub rep_thr: usize,

    #[arg(long = "threads")]
    pub threads: Option<usize>,
}

struct PlpCntRes {
    freq: Vec<f32>,
    depths: Vec<u32>,
    eq_cnts: Vec<u32>,
    seq: String,
}

impl PlpCntRes {
    pub fn new(freq: Vec<f32>, depths: Vec<u32>, eq_cnts: Vec<u32>, seq: String) -> Self {
        Self {
            freq,
            depths,
            eq_cnts,
            seq,
        }
    }
}

fn poly_ext(mut sbr_and_smc: SubreadsAndSmc, rep_thr: usize) -> SubreadsAndSmc {
    let smc_seq = &sbr_and_smc.smc.seq;
    let smc_seq_new = poly_ext_seq(smc_seq.as_bytes(), rep_thr);
    sbr_and_smc.smc.seq = smc_seq_new;
    sbr_and_smc
}

fn poly_ext_seq(seq: &[u8], rep_thr: usize) -> String {
    let mut pre_base = None;
    let mut cnt = 0;
    let mut res_seq = String::new();

    for &base in seq {
        if let Some(pre_base_) = pre_base {
            if pre_base_ == base {
                cnt += 1;
            } else {
                if cnt >= rep_thr {
                    res_seq.push(pre_base_ as char);
                }
                cnt = 1;
            }
        } else {
            cnt += 1;
        }
        res_seq.push(base as char);
        pre_base = Some(base);
    }

    if cnt >= rep_thr {
        res_seq.push(pre_base.unwrap() as char);
    }
    res_seq
}

fn polyn_ext_main(sbr_bam: &str, smc_bam: &str, threads: Option<usize>, rep_thr: usize) {
    tracing::info!("sorting sbr.bam {}", sbr_bam);
    let sorted_sbr = gskits::samtools::sort_by_tag(sbr_bam, "ch", None);

    tracing::info!("sorting smc.bam {}", smc_bam);
    let sorted_smc = gskits::samtools::sort_by_tag(smc_bam, "ch", None);

    let oup_filepath = format!("{}.polynext.bam", smc_bam.rsplit_once(".").unwrap().0);

    let (target2idx, idx2seq) = build_target2idx_and_idx2target(&sorted_smc);

    thread::scope(|s| {
        let sorted_sbr = &sorted_sbr;
        let sorted_smc = &sorted_smc;
        let target2idx = &target2idx;
        let idx2seq = &idx2seq;
        let (smc_sbr_sender, smc_sbr_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            asts::subreads_and_smc_generator(
                &sorted_sbr,
                &sorted_smc,
                &mm2::params::InputFilterParams::default(),
                smc_sbr_sender,
            )
        });

        let (ext_smc_sbr_sender, ext_smc_sbr_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            for sbr_and_smc in smc_sbr_recv {
                let ext_sbr_smc = poly_ext(sbr_and_smc, rep_thr);
                ext_smc_sbr_sender.send(ext_sbr_smc).unwrap();
            }
        });

        let threads = threads.unwrap_or(num_cpus::get());
        assert!(threads >= 4, "at lease 4 threads");
        let plp_therads = 5;

        let align_threads = threads - 2 - plp_therads;

        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..align_threads {
            let align_res_sender_ = align_res_sender.clone();
            let ext_smc_sbr_recv_ = ext_smc_sbr_recv.clone();
            s.spawn(move || {
                asts::align_worker(ext_smc_sbr_recv_, align_res_sender_, target2idx);
            });
        }
        drop(ext_smc_sbr_recv);
        drop(align_res_sender);

        let (plp_res_sender, plp_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..plp_therads {
            let align_res_recv_ = align_res_recv.clone();
            let plp_res_sender_ = plp_res_sender.clone();
            s.spawn(move || {
                for align_res in align_res_recv_ {
                    if align_res.records.len() == 0 {
                        continue;
                    }
                    let tid = align_res.records[0].tid() as usize;
                    let seq = idx2seq.get(&tid).unwrap().as_bytes();
                    let seq = poly_ext_seq(seq, rep_thr);

                    let plp_cnts = gskits::gsbam::plp_counts_from_records::plp_with_records_region(
                        &align_res.records,
                        Some(0),
                        Some(seq.len()),
                    );

                    // let plp_major_start = *plp_cnts.get_major().first().unwrap();
                    // let plp_major_end = *plp_cnts.get_major().last().unwrap() + 1;
                    // if plp_major_start != 0 || plp_major_end != seq.len() {
                    //     tracing::warn!(
                    //         "request plp_region:{}-{}, but got:{}-{}",
                    //         0,
                    //         seq.len(),
                    //         plp_major_start,
                    //         plp_major_end
                    //     );
                    // }

                    let cnts_array = Array2::from_shape_vec(
                        (10, plp_cnts.get_major().len()),
                        plp_cnts.get_cnts().to_vec(),
                    )
                    .unwrap();
                    let col_depth = cnts_array.sum_axis(Axis(0));
                    let mut major_depths = vec![];
                    let mut eq_cnts = vec![];
                    let mut freq = vec![];

                    let major_idx_in_plp_cnts = plp_cnts
                        .get_major()
                        .iter()
                        .zip(plp_cnts.get_minor().iter())
                        .enumerate()
                        .filter(|(_, (_, mi))| **mi == 0)
                        .map(|(idx, (&ma, _))| (ma, idx))
                        .collect::<HashMap<_, _>>();

                    (0..seq.len()).into_iter().for_each(|major_pos| {
                        let mut depth = 0;
                        let mut cnt = 0;
                        if let Some(idx) = major_idx_in_plp_cnts.get(&major_pos).map(|v| *v) {
                            depth = col_depth[idx];
                            cnt = cnts_array
                                    [[get_base_idx(seq.as_bytes()[major_pos], true), idx]]
                                    + cnts_array[[get_base_idx(seq.as_bytes()[major_pos], false), idx]];
                        }

                        major_depths.push(depth);
                        eq_cnts.push(cnt);
                        freq.push(if depth == 0 {0.0} else {cnt as f32 / depth as f32});

                    });

                    if seq.len() != freq.len() {
                        let major_points_set =
                            plp_cnts.get_major().iter().copied().collect::<HashSet<_>>();
                        let mut major_points = major_points_set.into_iter().collect::<Vec<_>>();
                        major_points.sort();

                        // let pos_set = plp_cnts
                        //     .get_major()
                        //     .iter()
                        //     .zip(plp_cnts.get_minor().iter())
                        //     .map(|(&ma, &mi)| (ma, mi))
                        //     .collect::<HashSet<_>>();

                        tracing::error!(
                            "seq_len:{}, qual_len:{}, major_pos:{:?}",
                            seq.len(),
                            freq.len(),
                            major_points
                        );
                        tracing::error!(
                            "major-minor:{:?}",
                            plp_cnts
                                .get_major()
                                .iter()
                                .zip(plp_cnts.get_minor().iter())
                                .collect::<Vec<_>>()
                        );

                        // for ma in plp_major_start..plp_major_end {
                        //     if !pos_set.contains(&(ma, 0)) {
                        //         tracing::error!("major init not found: {}", ma);
                        //     }
                        // }

                        // assert_eq!(major_points.len(), freq.len());

                        panic!();
                    }
                    plp_res_sender_
                        .send((tid, PlpCntRes::new(freq, major_depths, eq_cnts, seq)))
                        .unwrap();
                }
            });
        }
        drop(plp_res_sender);
        drop(align_res_recv);

        let pb = pbar::get_spin_pb("collecting plp info".to_string(), DEFAULT_INTERVAL);
        let plp_res = plp_res_recv
            .into_iter()
            .map(|plp_res| {
                pb.inc(1);
                plp_res
            })
            .collect::<HashMap<_, _>>();
        pb.finish();

        let mut target_bam_file = BamReader::from_path(smc_bam).unwrap();
        target_bam_file.set_threads(10).unwrap();
        let mut o_header = Header::from_template(target_bam_file.header());

        let mut hd = HeaderRecord::new(b"PG");
        hd.push_tag(b"ID", "polyNext")
            .push_tag(b"PN", "polyNext")
            .push_tag(b"CL", &command_line_str())
            .push_tag(b"VN", &env!("CARGO_PKG_VERSION"));

        if let Some(pp) = get_last_pg_from_bam_header(target_bam_file.header()) {
            hd.push_tag(b"PP", &pp);
        }
        o_header.push_record(&hd);

        let mut o_bam_file =
            BamWriter::from_path(&oup_filepath, &o_header, rust_htslib::bam::Format::Bam).unwrap();
        o_bam_file.set_threads(10).unwrap();

        let pb = pbar::get_spin_pb(
            format!("dump result to {}", oup_filepath),
            pbar::DEFAULT_INTERVAL,
        );

        for record in target_bam_file.records() {
            pb.inc(1);
            let mut record = record.unwrap();
            let record_ext = BamRecordExt::new(&record);
            let qname = record_ext.get_qname();
            let tid = target2idx.get(&qname).unwrap().0;

            if let Some(plp_cnts_res) = plp_res.get(&tid) {
                // let mut record_new = BamRecord::new();
                let qual = plp_cnts_res
                    .freq
                    .iter()
                    .map(|&f| quality_2_phreq(f, None))
                    .collect::<Vec<_>>();
                record.set(qname.as_bytes(), None, plp_cnts_res.seq.as_bytes(), &qual);
                record.remove_aux(b"rq").unwrap();
                record
                    .push_aux(
                        b"rq",
                        rust_htslib::bam::record::Aux::Float(gskits::phreq::phreq_list_2_quality(
                            &qual,
                        )),
                    )
                    .unwrap();

                record
                    .push_aux(b"dp", Aux::ArrayU32(AuxArray::from(&plp_cnts_res.depths)))
                    .unwrap();
                record
                    .push_aux(b"en", Aux::ArrayU32(AuxArray::from(&plp_cnts_res.eq_cnts)))
                    .unwrap();
            } else {
                tracing::warn!("no result for qname:{}", qname);
            }

            o_bam_file.write(&record).unwrap();
        }
        pb.finish();
    });
}

fn main() {
    let cli = Cli::parse();

    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero] [hour]:[minute]:[second]",
    )
    .unwrap();

    // let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);
    // let timer = tracing_subscriber::fmt::time::LocalTime::new(time_fmt);
    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);

    tracing_subscriber::fmt::fmt().with_timer(timer).init();

    polyn_ext_main(&cli.sbr, &cli.smc, cli.threads, cli.rep_thr);
}

#[cfg(test)]
mod test {
    use crate::poly_ext_seq;

    #[test]
    fn test_poly_ext_seq() {
        let seq = "AAACGA";
        let res_seq = poly_ext_seq(seq.as_bytes(), 3);
        assert_eq!(res_seq, "AAAACGA");

        let seq = "AAACGGA";
        let res_seq = poly_ext_seq(seq.as_bytes(), 2);
        assert_eq!(res_seq, "AAAACGGGA");

        let seq = "AAACGGAA";
        let res_seq = poly_ext_seq(seq.as_bytes(), 2);
        assert_eq!(res_seq, "AAAACGGGAAA");
    }
}
