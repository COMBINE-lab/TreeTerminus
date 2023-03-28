fn main() {
    println!("cargo:rerun-if-changed=phylip_src");
    cc::Build::new()
        .file("phylip_src/cons.c")
        .file("phylip_src/consense.c")
        .file("phylip_src/phylip.c")
        .include("phylip_src")
        .flag("-lm")
        .flag("-w")
        .compile("phylip_consense");
}
