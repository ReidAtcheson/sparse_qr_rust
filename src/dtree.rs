use crate::metis::MetisGraph;


pub struct DissectionTree{
     parents  : Vec<Option<usize>>,
     children : Vec<Option<(usize,usize)>>,
     levels : Vec<Vec<usize>>,
     nodes : Vec<Vec<usize>>
}



impl DissectionTree{
    ///Recursively splits input graph until every
    ///leaf has fewer than `maxnodes` nodes
    pub fn new(g : &MetisGraph,maxnodes : usize) -> Self{
        //First nodeset to start the tree off with
        let first : Vec<i64> = (0..g.nnodes()).map(|x|x as i64).collect();



        let (parents,nodes) = {
            let mut parents = Vec::<Option<usize>>::new();
            let mut nodes = Vec::<Vec<usize>>::new();
            let mut stack : Vec<(Vec<i64>,Option<usize>)> =  vec![(first,None)];
            //Populate `parents` and `nodes` 
            while let Some((ns,parent)) = stack.pop(){
                if ns.len()>maxnodes{
                    let subg=g.subgraph(&ns);
                    let (p1,p2,sep)=subg.split();

                    //The separator goes into the tree
                    nodes.push(sep.iter().map(|x|*x as usize).collect());
                    parents.push(parent);
                    //p1,p2 are partitions and go back into the stack,
                    //remembering the separator they split from.


                    //First we must recover p1,p2 in global
                    //graph vertex numbering.
                    let p1_global : Vec<i64> = p1.iter().map(|&id| ns[id as usize] ).collect();
                    let p2_global : Vec<i64> = p2.iter().map(|&id| ns[id as usize] ).collect();
                    let id=nodes.len()-1;
                    //Push nodes onto stack with their parents
                    stack.push( (p1_global,Some(id)) );
                    stack.push( (p2_global,Some(id)) );
                }
                else{
                    //This is the base case. Node sets now below the input
                    //threshold `maxnodes`
                    nodes.push(ns.iter().map(|&x|x as usize).collect());
                    parents.push(parent);
                }
            }
            (parents,nodes)
        };

        //Now populate `children`
        let children = {
            let mut children = vec![None;nodes.len()];
            for (i,sp) in parents.iter().enumerate(){
                if let Some(p) = sp{
                    children[*p]=match children[*p]{
                        None => Some((i,i)),
                        Some((k,_)) => {Some((k,i))}
                    }
                }
            }
            children
        };

        //Finally populate the `levels` array
        let levels = {
            let nlevels = {
                let mut nlevels=0;
                let mut stack = vec![(0,0)];
                while let Some((n,level)) = stack.pop(){
                    if let Some((c1,c2)) = children[n]{
                        stack.push((c1,level+1));
                        stack.push((c2,level+1));
                        nlevels = std::cmp::max(nlevels,level+1);
                    }
                }
                nlevels
            };

            let mut levels  : Vec<Vec<usize>> = vec![Vec::<usize>::new();nlevels];
            let mut stack = vec![(0,0)];
            levels[0].push(0);
            while let Some((n,level)) = stack.pop(){
                if let Some((c1,c2)) = children[n]{
                    stack.push((c1,level+1));
                    stack.push((c2,level+1));
                    levels[level].push(c1);
                    levels[level].push(c2);
                }
            }
            levels
        };

        DissectionTree { parents : parents, children : children,  levels : levels, nodes : nodes}
    }
}


#[cfg(test)]
mod tests {
    use crate::metis::MetisGraph;
    use crate::dtree::DissectionTree;


    #[test]
    fn dtree_construction() {
        let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
        let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
        let g = MetisGraph::new(xadj,adjncy); 
        let _dtree = DissectionTree::new(&g,3);
    }

    #[test]
    fn dtree_construction_stencil(){
        let mut xadj = Vec::<i64>::new();
        let mut adjncy = Vec::<i64>::new();
        let mut offs=0;

        let mx=64 as i64;
        let my=64 as i64;
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
        let _dtree = DissectionTree::new(&g,32);
    }

    #[test]
    fn dtree_construction_stencil_squared(){
        let mut xadj = Vec::<i64>::new();
        let mut adjncy = Vec::<i64>::new();
        let mut offs=0;

        let mx=64 as i64;
        let my=64 as i64;
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
        let _dtree = DissectionTree::new(&g2,32);
    }



}
