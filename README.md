# sparse_qr_rust
Simple nested-dissection based sparse QR factorization in Rust



# Dependencies

Uses [metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) to compute vertex separators. Not included in this repository.


I tried to set up an automated build step which grabs metis and its dependencies from github and then builds them with
the system `gcc`, pointing cargo to those libraries. If this step works then you don't need to do anything except
`cargo build`.


If it does not work then you need to build metis yourself and change `build.rs` to point to the resulting static
library `libmetis.a`.



# Example

You may find this file in `./examples/factorize_laplace2d.rs`

```rust
use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::numeric::SparseQR;
use sparse_qr_rust::gallery::laplace2d;


fn main(){
    let mx = 256;
    let my = 256;
    let m = mx*my;
    let maxnodes=200;
    type F=f64;
    //Assemble 2D Laplacian matrix with (mx*my) rows and (mx*my) columns
    let lap=laplace2d::<F>(mx,my);
    //Get squared graph
    let g = lap.to_metis_graph();
    let g2 = g.square();
    //Symbolic factorization
    let dtree = DissectionTree::new(&g2,maxnodes);
    //Numeric factorization, consumes input tree
    let mut fact = SparseQR::<F>::new(dtree,&lap);

    //Manufacture a solution to test the solver
    let soln : Vec<F> = (0..m).map(|x| (x as f64).sin()+2.0).collect();    
    let mut b = soln.clone();
    //Set b=A*soln where A is the Laplacian matrix
    lap.eval(&soln,&mut b);

    //Now go back and solve Ax=b and store result back into b
    fact.solve(&mut b);

    //Result `b` should be very close to `soln` by construction
    let mut errs : Vec<F> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).abs()/x.abs()).collect();
    errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
    let maxerr = errs.last().unwrap();
    println!("Maximum relative error: {}",maxerr);


}
```


# Next steps

Currently takes more memory than necessary because it conservatively represents every block in the nested dissection tree with dense matrices. This
makes using LAPACK easier but it's very wasteful. There are two likely steps to improve this:

1. Use "staircase" optimization for lower triangular parts 
    * these are the parts which get fed into xGEQRF
    * Stores the householder reflectors for the `Q` in `A=QR`
    * Also stores diagonal blocks for `R` which gets fed into xTRSM during solves
2. Use a hybrid sparse/dense format for upper triangular parts
    * This stores upper triangular part of `R`.
    * There is usually still a good deal of sparsity in these and only a few blocks truly benefit from dense representation
