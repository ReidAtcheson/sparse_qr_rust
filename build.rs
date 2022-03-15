use std::process::Command;


fn main() {
    let _output = Command::new("sh")
        .arg("./external/dl_install_metis.sh").output();
    println!("cargo:rustc-link-search=./external/lib")
}
