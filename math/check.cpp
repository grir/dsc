#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <set>
#include <vector>
#include <algorithm>
#include <immintrin.h>


#define i64 long long
using namespace std;


class Graph{
   public:
   enum RegularGraphs { R3_0 = 0, R3_1 = 1};
   i64 subs;
   int n;        // number of vertices
   int* data;    // adj matrix
   int* edgeMap; // edge numbers: 0,..., C(n,2)-1
   ///////////////////////////////////////////////
   Graph(){ 
      this->n = 0; 
      data = nullptr; 
      edgeMap = nullptr;
   }
   ///////////////////////////////////////////////
   void initByN(int n){
      this->n = n; 
      this->data = new int[n * n]; 
      this->edgeMap = new int[(n * (n - 1))/2]; 
      for(int i=0;i<n*n;i++)   
         data[i]=0;        
      int e = 0;
      for(int i=0;i<n;i++)   
         for(int j = i + 1;j < n;j++) {           
            edgeMap[ e ] = i * n + j;        
            e++;
         }   
   }
   ///////////////////////////////////////////////
   Graph(int n){ 
      initByN(n);
   }


   ///////////////////////////////////////////////
   Graph(mt19937_64& gen, int n, RegularGraphs cls){ // random regular (M<-2n-1)
      initByN(n);
      // K5:
      for(int i=0;i<5;i++)
        for(int j=i+1;j<5;j++)
          data[i * n + j] = data[j * n + i] = 1;
          
      // minus one edge
      data[4] = data[4 * n] = 0;
      
      // find v. with deg 3:
      //int sm = 2 * n; 
      int v0 = 5;

      while (v0 < n){

        int degs[n];
        
        for(int i=0;i<n;i++){
          int d = 0;
          for(int j=0;j<n;j++){
             d += data[i * n + j];
          }
          degs[i] = d;
                
        } 


        int v1,v2, v3, v4;
        for(int i=0; i<n;i++){
           for(int j=i+1; j<n;j++){
             if ((degs[i] + degs[j]) == 6) {
               v1 = i; v2 = j;           
             };           
             if ((degs[i] + degs[j]) == 8) {
               v3 = i; v4 = j;           
             };           
          }   
        }     
        
        data[v0 * n + v1] = data[v1 * n + v0] = 1;  
        data[v0 * n + v2] = data[v2 * n + v0] = 1;  
        data[v3 * n + v4] = data[v4 * n + v3] = 0;  
        data[v3 * n + v0] = data[v0 * n + v3] = 1;  
      
        v0++;
      }
          
   }


   ///////////////////////////////////////////////
   Graph(mt19937_64& gen, int n, int M){ // random GnM
      initByN(n);
      //cout << __LINE__ << endl;
      int Cn2 = n * (n - 1) / 2;
      int* temp = new int[Cn2];
     
      for(int i=0;i<Cn2;i++)
        temp[i] = i;

      //cout << __LINE__ << endl;

      uniform_int_distribution<int> rnde(0,Cn2-1);
      
      for(int i=0;i<M;i++){
        int rn = rnde(gen);
        int t = temp[i];
        temp[i] = temp[rn];
        temp[rn] = t;
        //cout << "-->" << temp[i] << endl;
      }
     
      for(int i = 0; i < M; i++)   
         data[ edgeMap[ temp[ i ] ] ] = 1;  
     
      //cout << __LINE__ << endl;
     
      for(int i=0;i<n;i++)   
         for(int j = i+1;j < n;j++) {           
            data[ j * n + i ] = data[ i * n + j ];            
         }       
         
      delete[] temp;   
    }  
   ///////////////////////////////////////////////
   Graph(string fname){
      ifstream fin(fname.c_str());
      int n;
      fin >> n;   
      initByN(n);

      for(int i=0;i<n*n;i++)
        fin >> data[i];

      fin.close();
   }
   void deepCopy(Graph& another){ //must be the same size!
      
      for(int i=0;i<n*n;i++)
        another.data[i] = this->data[i];
        
        
   }
   
   ///////////////////////////////////////////////
   
   void randomize(mt19937_64& gen, int numEdges){
       
     //unsigned long long  int  val;
     //cout << _rdrand64_step(&val) << endl; 
     int ns = n * n;
     uniform_int_distribution<int> rnde(0,ns - 1);
     for(int i=0;i<numEdges;i++){
         // choose r. edge
         int p1 = 0;         
         while (data[ p1 ] == 0)
            p1 = rnde(gen);
         // choose a place for new edge           
         int p2 = 0;         
         while ((data[ p2 ] != 0) || ((p2 % n) == (p2 / n)))
            p2 = rnde(gen);
         
         int i1 = p1 / n;
         int i2 = p2 / n;
         int j1 = p1 % n;
         int j2 = p2 % n;
         
         data[i1*n+j1] = data[j1*n+i1] = 0;       
         data[i2*n+j2] = data[j2*n+i2] = 1;       
     }
   
   }
   
   ///////////////////////////////////////////////
 
   friend ostream& operator<<(ostream& out, Graph& gr);
 
   ~Graph() { delete[] data; }
};

///////////////////////////////////////////////

ostream& operator<<(ostream& out, Graph& gr){
   int n = gr.n;
   for (int i=0;i<n;i++){
     for (int j=0;j<n;j++)
        out << gr.data[i*n+j] << " ";
     out << endl;    
   }
   return out;  
}


char getBit(i64 num, int bt){

    i64 sel = static_cast<i64>(1) << bt;
    if (num & sel) return 1;
    else return 0;

}
////////////////////////////////////////////
int numBits(i64 num, int n){ // returns sum of lowest n bits
    int sum = 0; 
    for(int i=0;i<n;i++)
       sum += getBit(num,i);
    return sum;

}

////////////////////////////////////////////////////////////////

int check4Subgraphs(Graph& gr, int minDeg){
   
   int n = gr.n;
   int* graph = gr.data; 
   i64 subs = 0;
   i64 nsubs = (static_cast<i64>(1) << n) - 1;
   //cout << "Subs " << nsubs << endl;
   //cout << "n " << n << endl;
   int minSubgr = n + 1;
   for (i64 ii = 1; ii < nsubs; ii++){
     int mins = n+1; // min v. degree in subgraph
     int numVerts = numBits(ii,n); // number of verts in the subgraph
     for (int i = 0; i < n; i++){
        if (!getBit(ii, i)) continue;
        int sumj = 0;
        for (int j = 0;j < n;j++){
           if (!getBit(ii, j)) continue;
           sumj += graph[i * n + j];
        }
        if (sumj < mins) {
          mins = sumj;            
          if (mins < minDeg) break;
        }        
     }
     if (mins >= minDeg) {
       subs = ii; 
       if(minSubgr > numVerts) 
          minSubgr = numVerts;
       //cout << subs << " verts -> " << numVerts << endl;
       //break;
     }  
   }  
   
   int lenSub = 0;
   for (int i=0;i<n;i++) 
     lenSub += getBit(subs, i);
   
   return minSubgr;   
   
}
const int MIL = 1<<20;

int bitSums[MIL];

void setBitSums(){
   for(int num = 0;num<MIL;num++){
      int nnum = num;
      int sum = 0;
      while (nnum>0){
        sum += nnum & 0x01;
        nnum >>= 1;
      }
      bitSums[num] = sum;
   }
}
////////////////////////////////////////////////////////////////
int getBitSum(i64 num){
   if (num < MIL) return bitSums[(int)num];
   int sum = 0;
   while (num > 0){
      sum += bitSums[num & (MIL-1)];
      num /= MIL;
   }
   return sum;
}

int check4SubgraphsF(Graph& gr, int minDeg, Graph& subgr){
   
   int n = gr.n;
   int* graphD = gr.data; 
   i64  graph[n];
   for(int i=0;i<n;i++){
      graph[i] = 0;
      i64 mul=1;
      for(int j=0;j<n;j++){
        graph[i] += mul * graphD[i*n+j];
        mul = mul << 1;
      }
   }
   //Graph subgr(n);  
   i64 subs = 0;
   i64 nsubs = (static_cast<i64>(1) << n) - 1;
   //cout << "Subs " << nsubs << endl;
   //cout << "n " << n << endl;
   int minSubgr = n + 1;
//   for (i64 ii = 1; ii < nsubs; ii++){
   for (i64 ii = nsubs; ii >= 0xF; ii--){

     int mins = n+1; // min v. degree in subgraph
     //int numVerts = numBits(ii,n); // number of verts in the subgraph
     int numVerts = getBitSum(ii & nsubs); // number of verts in the subgraph
     
     // next iteration if numVerts > size of current subgraph
     if ( (numVerts >= minSubgr) || (numVerts <= 3))  continue;
     //if (numVerts != (n-8))  continue;
     

     for (int i = 0; i < n; i++){
        if (!getBit(ii, i)) continue;
        int sumj = getBitSum(graph[i] & ii);
        //for (int j = 0;j < n;j++){
        //   if (!getBit(ii, j)) continue;
        //   sumj += graph[i * n + j];
        //}
        if (sumj < mins) {
          mins = sumj;            
          if (mins < minDeg) break;
        }        
     }
     if (mins >= minDeg) {
       gr.deepCopy(subgr);
       subgr.subs = ii; 
       if(minSubgr > numVerts) 
          minSubgr = numVerts;
       //cout << subs << " verts -> " << numVerts << endl;
       //break;
     } 
     //  if ((numVerts == n) 
   }  
   
   
   //cout << subgr << subs << endl;
   return minSubgr;   
   
}

///////////////////////////////////////////////////////

vector<int>& getMins(vector<int>& degs, 
                     vector<bool>& avbl, vector<int>& mins){
  int n = degs.size();
  // find min:
  int min = n+1;
  for(int i=0;i<n;i++)
     if ((min > degs[i]) && avbl[i]) min = degs[i];
  // find all mins:
  for(int i=0;i<n;i++)
     if ((min == degs[i]) && avbl[i]) mins.push_back(i);
     //if ( (min == (degs[i]-1) ) && avbl[i]) mins.push_back(i);
     
  return mins;
}

///////////////////////////////////////////////////////

void  remove(int* graph, int n, int v){
   for(int i = 0; i < n; i++){
     for(int j = 0;j < n; j++)
       if ((i == v) || (j == v)) 
         graph[ i * n + j ] = -1;
   }
}

///////////////////////////////////////////////////////

void  updateDegrees(int* graph, vector<int>& degs){
  int n = degs.size();
  for(int i = 0; i < n; i++){
     degs[ i ] = 0;
     for(int j = 0; j < n; j++)
        if (graph[i * n + j]>0) 
          degs[ i ] += graph[i * n + j];
   }
} 

///////////////////////////////////////////////////////
i64 toInt64(vector<bool>& v){
   //int n = v.size();	
   i64 res = 0;
   i64 twos = 1; 
   for(auto a : v){
      res = res  + (int) a * twos;
      twos <<= 1;
   }   
   return res;
}


///////////////////////////////////////////////////////


int heuristicCheck4Subgraphs(mt19937_64& gen, Graph& gr, int minDeg, Graph& subgr){
   
   int n = gr.n;
   int* graph = gr.data; 
   int minSubgr = n + 1;
   // (find degrees)   
   vector<int> degrees;
   vector<bool> avaible;
   int numAvaible = n;
   
   // init degrees and graph 
   for(int i = 0; i < n; i++){
     avaible.push_back( true );
     degrees.push_back( 0 );
     for(int j = 0; j < n; j++)
        degrees[ i ] += graph[i * n + j];
   }
   i64 subs = 0;
   while ( true ){
      vector<int> minList;
      getMins(degrees, avaible, minList);
      int sz = minList.size();
      int vidx = 0;
      if (sz>1){
         uniform_int_distribution<int> rnd(0,sz-1);    
         vidx = rnd(gen);
      }
         
      int v = minList[vidx];
      
      // update min subgraph:
      if ((degrees[ v ] >= minDeg) && (numAvaible < minSubgr)){
         minSubgr = numAvaible;
         gr.deepCopy(subgr);
         subgr.subs = toInt64(avaible);
      }   
      
      // remove from data the vertex v:
      remove(graph, n, v);
      avaible[v] = false;
      numAvaible--;
      if(numAvaible == 3) break; // the end of loop 
      
      //update degrees
      updateDegrees(graph, degrees); 
      
   }
   //cout << subgr << " " << subs << endl; 
   
   return minSubgr;   
   
}

///////////////////////////////////////////////////

i64 Cnk(int n, int k){
   i64 prod = 1;
   for(int i=0;i<k;i++)
      prod = (prod * (n-i)) /(i+1);
   return prod; 
}


////////////////////////////////////////////////////////////////

int main(){

  //Graph gr("gr10.dat");
  //cout << check4Subgraphs(gr, 3) << endl;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  
  mt19937_64 gen(seed);
  setBitSums();
  //cout << "Sums done " << endl;
  int n = 10;
  int M = 2 * n - 1;
  //int maxSub = 0;
  int opt = 0;
  int numIter = 100000;
  int devs = 0;
  /*for (int i=0;i<=32;i++ )
     cout << Cnk(32,i) << " ";*/
     
  //cout << "------------------" <<   endl;
  for (int i = 0; i < numIter; i++){

     //Graph gr(gen,n, Graph::RegularGraphs::R3_0);
     Graph gr(gen,n, M);

     //gr.randomize(gen, 20);
     Graph grx(n);
     gr.deepCopy(grx);

     // for subgraphs:
     Graph subgr(n);
     Graph subgrx(n);
     gr.deepCopy(subgr);
     gr.deepCopy(subgrx);
     
     
//      cout << gr << endl;  

     int csub = check4SubgraphsF(gr, 3, subgr);
    // int csubx = heuristicCheck4Subgraphs(gen, grx, 3, subgrx);
     
  //   if (csub > maxSub){
  //       maxSub = csub;
  
     if (csub == (n - 1))   
        cout << i << ", " << csub << endl 
        << hex << subgr.subs << endl << subgr;  
     opt++;
 //    }
   /*  else {
       cout << "------------------" <<   endl;
       cout << "i :" << i << " Exact --> " << subgr << hex << subgr.subs << " " << dec << csub << endl;
       cout << "heuristic --> " << subgrx << hex << subgrx.subs << " " << dec << csubx << endl;
    }
       */
//       devs += (csubx - csub);
       
//       if (csub > csubx)   cout << "Something wrong" << endl;  
       
  //   }    
  //  if (i % 500 == 0) cout << " --> " << i << " " << opt << " devs: " << devs*1.0/(i+1) << endl;

  }
  
}
