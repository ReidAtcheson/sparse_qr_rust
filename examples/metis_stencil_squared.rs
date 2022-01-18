///Computes simple stencil graph, squares it,
///and then prints out graphviz dot file
///for both original and squared

use sparse_qr_rust::metis::MetisGraph;

fn main() {

    let mut xadj = Vec::<i64>::new();
    let mut adjncy = Vec::<i64>::new();
    let mut offs=0;

    let mx=8 as i64;
    let my=8 as i64;
    let id = |ix : i64,iy : i64|{iy + my*ix};
    xadj.push(offs);
    for ix in 0..mx{
        for iy in 0..my{

            if ix>0{
                adjncy.push(id(ix-1,iy));
                offs+=1;
            }
            if ix<mx-1{
                adjncy.push(id(ix+1,iy));
                offs+=1;
            }
            if iy>0{
                adjncy.push(id(ix,iy-1));
                offs+=1;
            }
            if iy<my-1{
                adjncy.push(id(ix,iy+1));
                offs+=1;
            }
            xadj.push(offs);
        }
    }
    let g = MetisGraph::new(xadj,adjncy);
    let g2 = g.square();
    g.graphviz("stencil.dot".to_string()).unwrap();
    g2.graphviz("stencil2.dot".to_string()).unwrap();
}
