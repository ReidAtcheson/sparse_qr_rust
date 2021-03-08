use sparse_qr_rust::sparse::CSCSparse;
use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::metis::MetisGraph;
use sparse_qr_rust::numeric::SparseQR;
use sparse_qr_rust::gallery::laplace2d;



fn main(){
    let mx = 512;
    let my = 512;
    let maxnodes=500;
    type F=f64;
    //Assemble matrix
    let lap=laplace2d::<F>(mx,my);
    //Get squared graph
    let g = lap.to_metis_graph();
    let g2 = g.square();
    //Symbolic factorization
    let dtree = DissectionTree::new(&g2,maxnodes);
    //Numeric factorization, consumes input tree
    let _fact = SparseQR::<F>::new(dtree,&lap);
}
