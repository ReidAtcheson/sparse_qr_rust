use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::string::String;
use std::fs::File;
use std::io::Write;

#[link(name = "metis", kind = "static")]
extern{
    pub fn METIS_ComputeVertexSeparator(nvtxs : *const i64,xadj : *const i64,adjncy : *const i64,
                                    vwgt : *const i64, options : *const i64, sepsize : *mut i64,part : *mut i64) -> i32;
}



pub struct MetisGraph{
    xadj : Vec<i64>,
    adjncy : Vec<i64> 
}


impl MetisGraph{
    pub fn nnodes(&self) -> usize { self.xadj.len()-1 }
    pub fn panic_if_invalid(&self) -> (){
        assert!(self.xadj.len()>0);
        assert!(self.adjncy.len()>0);
        assert!(*self.xadj.last().unwrap()==((self.adjncy.len()) as i64));
        for edge in self.adjncy.iter(){
            assert!(*edge<(self.xadj.len()-1) as i64);
        }
        //TODO: Replace with "is_sorted" when it merges
        for i in 1..self.xadj.len(){
            assert!(self.xadj[i-1]<=self.xadj[i]);
        }

        //Metis doesn't like self-connections, make sure there are none
        for i in 1..self.xadj.len(){
            let beg=self.xadj[i-1];
            let end=self.xadj[i];
            for e in self.adjncy[beg as usize..end as usize].iter().cloned(){
                assert!( (i-1) != e as usize);
            }
        }

        //Make sure graph is structurally symmetric
        for i in 1..self.xadj.len(){
            let beg=self.xadj[i-1];
            let end=self.xadj[i];
            for e in self.adjncy[beg as usize..end as usize].iter().cloned(){
                let beg2=self.xadj[e as usize];
                let end2=self.xadj[(e+1) as usize];
                let n = (i-1) as i64;
                assert!(self.adjncy[beg2 as usize..end2 as usize].contains(&n));
            }
        }
    }

    pub fn new(xadj : Vec<i64>, adjncy : Vec<i64>) -> Self{
        let out = MetisGraph { xadj : xadj, adjncy : adjncy };
        out.panic_if_invalid();
        out
    }

    ///Computes subgraph by taking a list of
    ///vertices in the original global numbering
    ///and producing graph where an edge exists between
    ///the new vertices if the edge existed between them
    ///in the original graph.
    pub fn subgraph(&self,pids : &Vec<i64>) -> Self{
        //These will hold data for the subgraph
        let mut xadj = Vec::<i64>::new();
        xadj.push(0);
        let mut adjncy = Vec::<i64>::new();

        //This inverts the relationship between local ids -> global ids
        let mut idmap  = BTreeMap::<i64,i64>::new();
        for (j,p) in pids.iter().enumerate(){
            idmap.insert(*p,j as i64);
        }


        let mut offs=0 as i64;
        for p in pids.iter(){
            //Get all the edges connected to vertex `p`
            let beg = self.xadj[*p as usize];
            let end = self.xadj[(*p+1) as usize];
            //Loop over all these edges, determine if 
            //the edge connects to a vertex in the input set.
            //If so, it goes in the subgraph in the new local numbering.
            for j in (beg as usize)..(end as usize){
                let gedge=&self.adjncy[j];
                if let Some(ledge) = idmap.get(gedge){
                    adjncy.push(*ledge); 
                    offs+=1;
                }
            }
            xadj.push(offs);
        }
        MetisGraph::new(xadj,adjncy)
    }
    pub fn split(&self) -> (Vec<i64>,Vec<i64>,Vec<i64>){
        const PART0 : i64 = 0;
        const PART1 : i64 = 1;
        const SEPARATOR : i64 = 2;
        const METIS_OK : i32 = 1;
        const METIS_ERROR_INPUT : i32 = -2;
        const METIS_ERROR_MEMORY : i32 = -3;
        const METIS_ERROR : i32 = -4;
        

        let part  = {
            let nvtxs : i64 = (self.xadj.len()-1) as i64;
            let mut part : Vec<i64> = vec![0;nvtxs as usize];
            let mut sepsize : i64 = 0;
            let info : i32 = 
            unsafe{ 
                METIS_ComputeVertexSeparator(&nvtxs,self.xadj.as_ptr(),self.adjncy.as_ptr(),std::ptr::null(),std::ptr::null(),&mut sepsize,part.as_mut_ptr())
            };
            assert!(sepsize != 0);
            assert!(info != METIS_ERROR_INPUT);
            assert!(info != METIS_ERROR_MEMORY);
            assert!(info != METIS_ERROR);
            assert_eq!(info,METIS_OK);
            part
        };

        let p0 : Vec<i64> = part.iter().enumerate().filter(|&(_i,v)| *v == PART0).map(|(i,_v)|i as i64).collect();
        let p1 : Vec<i64> = part.iter().enumerate().filter(|&(_i,v)| *v == PART1).map(|(i,_v)|i as i64).collect();
        let psep : Vec<i64> = part.iter().enumerate().filter(|&(_i,v)| *v == SEPARATOR).map(|(i,_v)|i as i64).collect();

        //These cases can actually happen in practice
        //but handling them adds a lot of complications
        //so for now I panic when they happen.
        assert!(p0.len()>0);
        assert!(p1.len()>0);
        assert!(psep.len()>0);
        assert!(p0.len()+p1.len()+psep.len()==self.xadj.len()-1);

        (p0,p1,psep)
    }

    ///Compute the squared graph
    pub fn square(&self) -> Self{
        let mut offs=0;
        let mut xadj = Vec::<i64>::new();
        xadj.push(offs);
        let mut adjncy = Vec::<i64>::new();
        for i in 1..self.xadj.len(){
            let beg=self.xadj[i-1];
            let end=self.xadj[i];
            let mut s : BTreeSet<i64> = self.adjncy[beg as usize..end as usize].iter().cloned().collect();
            for e1 in beg..end{
                let e = self.adjncy[e1 as usize];
                let beg2 = self.xadj[e as usize];
                let end2 = self.xadj[e as usize + 1 as usize];
                let t : BTreeSet<i64> = self.adjncy[beg2 as usize..end2 as usize].iter().cloned().collect();
                s=s.union(&t).cloned().collect();
                //Avoid self connections because they make METIS unhappy
                s.remove(&((i-1) as i64));
            }
            offs+=s.len() as i64;
            xadj.push(offs);
            for x in s.iter(){
                adjncy.push(*x);
            }
        }
        MetisGraph::new(xadj,adjncy)
    }

    pub fn graphviz(&self,fname : String) -> std::io::Result<()>{
        let mut file = File::create(fname)?;
        write!(&mut file,"graph g {{ \n")?;
        for i in 1..self.xadj.len(){
            let beg=self.xadj[i-1];
            let end=self.xadj[i];
            for e in self.adjncy[beg as usize .. end as usize].iter(){
                if i-1<=*e as usize{
                    write!(&mut file,"  {} -- {}\n",i-1,e)?;
                }

            }
        }
        write!(&mut file,"}} \n")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::metis::MetisGraph;
    use std::collections::BTreeSet;
    #[test]
    fn metis_graph_construction() {
        let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
        let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
        let _g = MetisGraph::new(xadj,adjncy); 
    }

    #[test]
    fn metis_graph_split(){
        let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
        let nnodes=xadj.len()-1;
        let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
        let g = MetisGraph::new(xadj,adjncy); 
        let (p1,p2,sep) = g.split();

        let s : BTreeSet<i64> = (0..nnodes as i64).collect();
        let s1 : BTreeSet<i64> = p1.iter().map(|&x|x).collect();
        let s2 : BTreeSet<i64> = p2.iter().map(|&x|x).collect();
        let s3 : BTreeSet<i64> = sep.iter().map(|&x|x).collect();

        //Sets s1,s2,s3 should be disjoint from each other
        assert!(s1.is_disjoint(&s2));
        assert!(s1.is_disjoint(&s3));
        assert!(s2.is_disjoint(&s3));
        let s12 : BTreeSet<i64> = s1.union(&s2).map(|&x|x).collect();
        let s123 : BTreeSet<i64> = s3.union(&s12).map(|&x|x).collect();
        assert!(s.is_subset(&s123));
        assert!(s123.is_subset(&s));


    }

    #[test]
    fn metis_graph_subgraph(){
        let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
        let _nnodes=xadj.len()-1;
        let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
        let g = MetisGraph::new(xadj,adjncy); 
        let (p1,p2,_sep) = g.split();
        let _g1 = g.subgraph(&p1);
        let _g2 = g.subgraph(&p2);
    }

    #[test]
    fn metis_graph_squared(){
        let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
        let _nnodes=xadj.len()-1;
        let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
        let g = MetisGraph::new(xadj,adjncy); 
        let _g2 = g.square();
    }


}
