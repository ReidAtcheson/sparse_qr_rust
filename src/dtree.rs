use std::fs::File;
use std::io::Write;
use std::string::String;
use crate::metis::MetisGraph;
use std::collections::BTreeSet;


pub struct DissectionTree{
     pub parents  : Vec<Option<usize>>,
     pub children : Vec<Option<(usize,usize)>>,
     pub levels : Vec<Vec<usize>>,
     pub nodes : Vec<Vec<usize>>
}



impl DissectionTree{
    pub fn panic_if_invalid(&self) -> () {
        assert_eq!(self.levels[0].len(),1);
        assert_eq!(self.parents[0],None);
        assert_eq!(self.parents.len(),self.children.len());
        assert_eq!(self.nodes.len(),self.parents.len());

        for (j,c) in self.children.iter().cloned().enumerate(){
            if let Some((c1,c2)) = c{
                let pc1 = self.parents[c1].unwrap();
                let pc2 = self.parents[c2].unwrap();
                assert_eq!(pc1,pc2);
                assert_eq!(pc1,j);
                assert_eq!(pc2,j);
            }
        }
        for (j,pm) in self.parents.iter().cloned().enumerate(){
            if let Some(p) = pm{
                let (c1,c2) = self.children[p].unwrap();
                assert!(c1 == j || c2 == j);
            }
        }

        for j in 0..self.levels.len()-1{
            for n in self.levels[j].iter().cloned(){
                let next_level : BTreeSet<usize> = self.levels[j+1].iter().cloned().collect();
                if let Some((c1,c2)) = self.children[n]{
                    assert!(next_level.contains(&c1));
                    assert!(next_level.contains(&c2));
                }
            }
        }
    }
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
                    if sep.len() == 0{
                        nodes.push(ns.iter().map(|&x|x as usize).collect());
                        parents.push(parent);
                        continue;
                    }


                    //The separator goes into the tree
                    nodes.push(sep.iter().map(|&id|ns[id as usize] as usize).collect());
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
            let mut children : Vec<Option<(usize,usize)>> = vec![None;nodes.len()];
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
                nlevels+1
            };

            let mut levels  : Vec<Vec<usize>> = vec![Vec::<usize>::new();nlevels];
            let mut stack = vec![(0,0)];
            levels[0].push(0);
            while let Some((n,level)) = stack.pop(){
                if let Some((c1,c2)) = children[n]{
                    stack.push((c1,level+1));
                    stack.push((c2,level+1));
                    levels[level+1].push(c1);
                    levels[level+1].push(c2);
                }
            }
            levels
        };

        let dtree = DissectionTree { parents : parents, children : children,  levels : levels, nodes : nodes};
        dtree.panic_if_invalid();
        dtree
    }

    pub fn get_permutation(&self) -> Vec<usize>{
        let nnodes=self.nodes.iter().map(|n|n.len()).fold(0,|acc,x|acc+x);
        let mut p = vec![0 as usize;nnodes];
        let mut it=0;
        for level in self.levels.iter().rev(){
            for n in level.iter(){
                for i in self.nodes[*n].iter(){
                    p[it]=*i;
                    it+=1;
                }
            }
        }
        p
    }
    //Plot tree to dotfile using parent relationships
    pub fn to_graphviz_parents(&self,fname : String) -> std::io::Result<()>{
        let mut file = File::create(fname)?;
        write!(&mut file,"digraph g {{ \n")?;

        for (i,sp) in self.parents.iter().enumerate(){
            if let Some(p) = sp{
                write!(&mut file,"  {} -> {}\n",i,p)?;
            }
        }
        write!(&mut file,"}} \n")?;
        Ok(())
    }

    //Plot tree to dotfile using children relationships
    pub fn to_graphviz_children(&self,fname : String) -> std::io::Result<()>{
        let mut file = File::create(fname)?;
        write!(&mut file,"digraph g {{ \n")?;

        for (i,c) in self.children.iter().enumerate(){
            if let Some((c1,c2)) = c{
                write!(&mut file,"  {} -> {}\n",i,c1)?;
                write!(&mut file,"  {} -> {}\n",i,c2)?;
            }
        }
        write!(&mut file,"}} \n")?;
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;
    use crate::metis::MetisGraph;
    use crate::dtree::DissectionTree;
    use crate::gallery::laplace2d;

    //Node sets should be a partition
    #[test]
    fn dtree_nodes_partition(){
        let mx=32;
        let my=32;
        let m=mx*my;
        let a = laplace2d::<f64>(mx,my);
        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,100);

        let s1 : BTreeSet<usize> = (0..m).collect();
        let mut s2 = BTreeSet::<usize>::new();
        for n in dtree.nodes.iter(){
            for i in n.iter(){
                s2.insert(*i);
            }
        }
        let tm = dtree.nodes.iter().map(|n|n.len()).fold(0,|acc,x|acc+x);

        assert_eq!(m,tm);
        assert_eq!(s2.len(),m);
        assert_eq!(s1,s2);
    }



    //Nested dissection tree should produce a permutation
    #[test]
    fn dtree_permutation() {
        let mx=32;
        let my=32;
        let m=mx*my;
        let a = laplace2d::<f64>(mx,my);
        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,100);
        let perm = dtree.get_permutation();        
        //Length of permutation array should be equal to number of equations (rows)
        //in the input sparse matrix.
        assert_eq!(perm.len(),m);
        let s1 : BTreeSet<usize> = (0..m).collect();
        let s2 : BTreeSet<usize> = perm.iter().cloned().collect();
        //Check for no duplicate entries in the permutation array
        assert_eq!(s2.len(),perm.len());
        //Finally, all indices should be represented in the permutation array
        assert_eq!(s1,s2);
    }


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
        let _dtree = DissectionTree::new(&g,64);
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
