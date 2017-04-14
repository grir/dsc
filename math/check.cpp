#include <fstream>
#include <iostream>

#define i64 long long
using namespace std;


class Graph{
   public:
   
   int n;
   
   int* data;
  
   int* edgeMap; 
   
   Graph(){ this->n = 0; data = nullptr; edgeMap = nullptr;}
   
   initByN(int n){
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
   
   Graph(int n){ 
      initByN(int n);
   }
   
   Graph(string fname){
      ifstream fin(fname.c_str());
      int n;
      fin >> n;   
      initByN(n);

      for(int i=0;i<n*n;i++)
        fin >> data[i];

      fin.close();
   }
   
   int edgIdx(int edge){
      
   
   
   
   
   }
   
   ~Graph() {delete[] data;}
};

////////////////////////////////////////////////////////////////

char getBit(i64 num, int bt){

    i64 sel = static_cast<i64>(1) << bt;
    if (num & sel) return 1;
    else return 0;

}

////////////////////////////////////////////////////////////////

int check4Subgraphs(Graph& gr, int minDeg){
   
   int n = gr.n;
   int* graph = gr.data; 
   i64 subs = 0;
   i64 nsubs = (static_cast<i64>(1) << n);//-1;
   cout << "Subs " << nsubs << endl;
   cout << "n " << n << endl;

   for (i64 ii = 1; ii < nsubs; ii++){
     int mins = n+1;
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
     if (mins>=minDeg) {
       subs = ii; 
       cout << subs << endl;
       break;
     }  
   }  
   
   int lenSub = 0;
   for (int i=0;i<n;i++) 
     lenSub += getBit(subs, i);
   
   return lenSub;   
   
}

////////////////////////////////////////////////////////////////

int main(){

  Graph gr("g11.dat");
    
  cout << "lenSub " << check4Subgraphs(gr, 3) << endl;  

}

