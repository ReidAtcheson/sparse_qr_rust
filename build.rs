use std::process::Command;
use std::fs;
use std::path::PathBuf;
use std::env;

fn main() {
    let args : Vec<String> = env::args().collect();
    let relfile = PathBuf::from(args[0].clone());
    let absfile = fs::canonicalize(&relfile);
    if let Ok(afile) = absfile{
        if let Some(astr) = afile.to_str(){
            let abspath = (&astr[0..astr.len()-5]).to_string();
            let build = abspath.clone() + &"/external/dl_install_metis.sh".to_string();
            let lib = abspath.clone() + &"/external/lib".to_string();

            let output = Command::new("sh").arg(build).output().unwrap();
            //println!("{:?}",output);
            println!("cargo:rustc-link-search={}",lib);
            println!("cargo:rustc-link-lib=blas");
            println!("cargo:rustc-link-lib=lapack");
        }
    }
}
