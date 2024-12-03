use std::collections::HashMap;

use rust_htslib::bam::Read;

type IdxQlen = (usize, usize);

pub fn build_target2idx_and_idx2target(smc_bam: &str) -> (HashMap<String, IdxQlen>, HashMap<usize, String>) {
    let mut reader = rust_htslib::bam::Reader::from_path(smc_bam).unwrap();
    reader.set_threads(10).unwrap();

    let mut target2idx = HashMap::new();
    let mut idx2seq = HashMap::new();

    for (idx, record) in reader.records().enumerate() {
        let record = record.unwrap();

        let qname = unsafe { String::from_utf8_unchecked(record.qname().to_owned()) };

        let qlen = record.seq_len();
        target2idx.insert(qname, (idx, qlen));
        idx2seq.insert(idx, unsafe{String::from_utf8_unchecked(record.seq().as_bytes())});
    }
    (target2idx, idx2seq)
}
