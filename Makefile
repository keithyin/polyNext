build:
	cargo build --release

install:
	cp target/release/polyNext /usr/bin/

clean:
	rm -rf target