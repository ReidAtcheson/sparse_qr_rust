use crate::metis::MetisGraph;


pub struct DissectionTree{
     parents  : Vec<Option<usize>>,
     children : Vec<Option<(usize,usize)>>,
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

        DissectionTree { parents : parents, children : children,  nodes : nodes}
    }
}
