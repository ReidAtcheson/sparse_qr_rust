use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::numeric::SparseQR;
use sparse_qr_rust::gallery::laplace1d;



fn main(){
    let m = 10000;
    let maxnodes=200;
    type F=f64;
    //Assemble matrix
    let lap=laplace1d::<F>(m);
    //Get squared graph
    let g = lap.to_metis_graph();
    let g2 = g.square();
    //Symbolic factorization
    let dtree = DissectionTree::new(&g2,maxnodes);
    //Numeric factorization, consumes input tree
    let _fact = SparseQR::<F>::new(dtree,&lap);
}
