# sparse_qr_rust
Simple nested-dissection based sparse QR factorization in Rust



# Dependencies

Uses [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) to compute vertex separators. Not included in this repository.
To use, do the following.
  1. Download metis follow install instructions.
    * Be sure you switch datatypes to 64 bits as they are 32 bits by default. Install instructions tell you how to do this.
  2. Change `build.rs` to point to the directory that you intalled metis
  3. That should be all you need to do
