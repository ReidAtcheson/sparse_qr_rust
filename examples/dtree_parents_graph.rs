use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::gallery::laplace2d;



fn main(){
    let mx = 64;
    let my = 64;
    let maxnodes=200;
    type F=f64;
    //Assemble matrix
    let lap=laplace2d::<F>(mx,my);
    //Get squared graph
    let g = lap.to_metis_graph();
    let g2 = g.square();
    //Symbolic factorization
    let dtree = DissectionTree::new(&g2,maxnodes);
    dtree.to_graphviz_parents("dtree_parents.dot".to_string()).ok();
}
