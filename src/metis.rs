use std::collections::BTreeMap;

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
        assert!(*self.xadj.last().unwrap()==((self.adjncy.len()+1) as i64));
        for edge in self.adjncy.iter(){
            assert!(*edge<(self.xadj.len()-1) as i64);
        }
        //TODO: Replace with "is_sorted" when it merges
        for i in 1..self.xadj.len(){
            assert!(self.xadj[i-1]<=self.xadj[i]);
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
                match idmap.get(&self.adjncy[j]){
                    Some(q) => {adjncy.push(*q); offs+=1;}
                    None => ()
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

        let part  = {
            let nvtxs : i64 = (self.xadj.len()-1) as i64;
            let mut part : Vec<i64> = vec![0;nvtxs as usize];
            let mut sepsize : i64 = 0;
            unsafe{ 
                METIS_ComputeVertexSeparator(&nvtxs,self.xadj.as_ptr(),self.adjncy.as_ptr(),std::ptr::null(),std::ptr::null(),&mut sepsize,part.as_mut_ptr());
            }
            assert!(sepsize==3);
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
}


