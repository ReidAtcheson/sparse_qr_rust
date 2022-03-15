use std::process::Command;
use std::path::Path;
use std::env;

fn main() {
    let dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let build = Path::new(&dir).join("external").join("dl_install_metis.sh");
    let lib = Path::new(&dir).join("external").join("lib");
    let _output = Command::new("sh").arg(build.to_str().unwrap()).output().unwrap();
    println!("cargo:rustc-link-search={}",lib.display());
    println!("cargo:rustc-link-lib=blas");
    println!("cargo:rustc-link-lib=lapack");
}
