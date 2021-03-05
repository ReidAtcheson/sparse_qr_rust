use crate::dtree::DissectionTree;
use crate::sparse::CSCSparse;
use crate::utility::Lapack;
use std::collections::BTreeMap;
use num_traits::Num;


pub struct SparseQR<F>{
    //Tree data
    parents  : Vec<Option<usize>>,
    children : Vec<Option<(usize,usize)>>,
    levels : Vec<Vec<usize>>,
    nodes : Vec<Vec<usize>>,
    //Factorization data
    nrows : usize,
    ncols : usize,

    //Upper triangular blocks for a block of columns.
    //Blocks of columns corresponding to leafs in
    //the nested dissection tree have no upper
    //triangular blocks and so may have zero
    //length
    triu : Vec<Vec<usize>>,
    //Lower triangular blocks for a block of columns.
    //Including diagonal block.
    tril : Vec<Vec<usize>>,

    //Temporary storage for submatrices.
    //Option type here used to deallocate
    //submatrices when no longer needed.
    temp : Vec<BTreeMap<usize,Vec<F>>>,


    //Packed storage so that solve phase can proceed
    //with only two calls to LAPACK.
    //Upper triangular block. 
    triu_num : Vec<Vec<F>>,
    //Lower triangular block.
    tril_num : Vec<Vec<F>>
} 


impl <F : Lapack+Num+Copy> SparseQR<F>{
    pub fn new(dtree : DissectionTree, mat : &CSCSparse<F>) -> Self{

        let nnodes=dtree.nodes.len();
        //Propagate lower triangular parts down tree
        let mut tril : Vec::<Vec<usize>> = vec![Vec::<usize>::new();nnodes];
        for node in 0..nnodes{
            //Lower triangular part contains the diagonal block as well.
            tril[node].push(node);
            //Now propagate down tree
            let mut stack : Vec<usize> = vec![node];
            while let Some(p) = stack.pop(){
                if let Some((c1,c2)) = dtree.children[p]{
                    tril[c1].push(node);
                    tril[c2].push(node);
                    stack.push(c1);
                    stack.push(c2);
                }
            }
        }
        //Propagate upper triangular parts 
        let mut triu : Vec::<Vec<usize>> = vec![Vec::<usize>::new();nnodes];
        for node in 0..nnodes{
            let mut maybep = dtree.parents[node];
            //Loop through the path to the top node
            while let Some(p) = maybep{
                triu[p].push(node);
                maybep=dtree.parents[p];
            }
        }
        //Sort the upper and lower triangular block lists
        //This step isn't strictly necessary but it does
        //make some logic cleaner later.
        //
        //If we sort these lists then we get the following
        //properties:
        //
        //tril[0] is always a diagonal block
        //tril[i+1] is always a block corresponding to one level down from tril[i]
        //triu[i-1] is always a block corresponding to _at most_ one level up from triu[i]
        for mut tri in tril.iter_mut(){
            tri.sort();
        }
        for mut tri in triu.iter_mut(){
            tri.sort();
        }

        //Assemble initial matrices into the 
        //temporary matrix list
        let mut temp : Vec<BTreeMap<usize,Vec<F>>> = vec![BTreeMap::<usize,Vec<F>>::new();nnodes];
        for node in 0..nnodes{
            let cols  = &dtree.nodes[node];
            //Assemble lower triangular submatrices
            for tri in tril[node].iter(){
                let rows = &dtree.nodes[*tri];
                let a = mat.slice_copy(rows,cols);
                temp[node].insert(*tri,a);
            }
            //Assemble upper triangular submatrices
            for tri in triu[node].iter(){
                let rows = &dtree.nodes[*tri];
                let a = mat.slice_copy(rows,cols);
                temp[node].insert(*tri,a);
            }
        }
        //Begin factorization.


        //Loop from bottom level of dissection tree to top level
        let mut tril_num : Vec::<Vec<F>> = vec![Vec::<F>::new();nnodes];
        for level in dtree.levels.iter().rev(){
            for node in level.iter().cloned(){
                //Gather lower triangular blocks into a stacked matrix
                let ncols=dtree.nodes[node].len();
                let nrows=tril[node].iter().map(|&k|dtree.nodes[k].len()).fold(0,|acc,x|acc+x);
                tril_num[node]=vec![F::zero();nrows*ncols];
                let trilmat=&mut tril_num[node];
                let nmap=&temp[node];
                let mut offs=0;
                for tri in tril[node].iter().cloned(){
                    let m = nmap.get(&tri).unwrap();
                    let slice=&mut trilmat[offs..offs+m.len()];
                    for (x,y) in slice.iter_mut().zip(m.iter()){
                        *x=*y;
                    }
                    offs+=m.len();
                }
                //TODO: query optimal workspace size. add tau vectors to the numeric qr struct
                let mut tau : Vec<F> = vec![F::zero();ncols];
                let nb=32;
                F::xgeqrf(nrows as i32,ncols as i32,&mut trilmat,nrows as i32,&tau,
            }
        }







        let triu_num = Vec::<Vec<F>>::new();
        SparseQR { parents : dtree.parents, children : dtree.children, levels : dtree.levels, nodes : dtree.nodes,
        nrows : mat.get_nrows(), ncols : mat.get_ncols(), triu : triu, tril : tril, temp : temp, triu_num : triu_num, tril_num : tril_num}
    }
}



#[cfg(test)]
mod tests {
    use crate::sparse::CSCSparse;
    use crate::metis::MetisGraph;
    use crate::dtree::DissectionTree;
    use crate::numeric::SparseQR;


    #[test]
    fn construct_stencil_factorization(){
        let mut offsets = Vec::<usize>::new();
        let mut rids = Vec::<usize>::new();
        let mut vals = Vec::<f64>::new();
        let mut offs=0;

        let mx=32 as usize;
        let my=32 as usize;
        let nrows=mx*my;
        let ncols=mx*my;
        let id = |ix : usize,iy : usize|{iy + my*ix};
        offsets.push(offs);
        for ix in 0..mx{
            for iy in 0..my{
                rids.push(id(ix,iy));
                vals.push(1.0);
                offs+=1;
                if ix>0{
                    rids.push(id(ix-1,iy));
                    vals.push(1.0);
                    offs+=1;
                }
                if ix<mx-1{
                    rids.push(id(ix+1,iy));
                    vals.push(1.0);
                    offs+=1;
                }
                if iy>0{
                    rids.push(id(ix,iy-1));
                    vals.push(1.0);
                    offs+=1;
                }
                if iy<my-1{
                    rids.push(id(ix,iy+1));
                    vals.push(1.0);
                    offs+=1;
                }
                offsets.push(offs);
            }
        }

        let a = CSCSparse::new(nrows,ncols,offsets,rids,vals);
        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,30);
        let fact = SparseQR::new(dtree,&a);
    }


}
