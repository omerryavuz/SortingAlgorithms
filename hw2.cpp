//============================================================================
// Name        : hw2.cpp
// Author      : Omer Yavuz
// Version     : ceng315-hw2
//============================================================================





/* I use two bucketsort function : 
 * bucketsort is used to sort speeds (0.xxx)
 * bucketsort2 is used to sort points (integers) 
 * 
 * Time and speed are stored as double and sorting functions work without rounding. 
 * 
 * 
 */


#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <string.h>
#include <math.h>

using namespace std;

#define N 5

int numberOfAnt;

typedef pair<double, int> Ant_pair;

class Ant{
	public:
		int id;
		int point;
		vector<double> speed ;
		vector<double> time ;
		Ant(){
			point = 0;
			id = 0;
		}
};


vector<double> splitString(string line){
	vector<double> list;
	string temp = "";
	for(unsigned int i=0 ; i<line.length() ; i++){
		if(line[i] == ' '){
			double valueOfstr = atof(temp.c_str());
			list.push_back(valueOfstr);
			temp = "";
		}
		else{
			temp += line[i];
		}
	}
	list.push_back(atof(temp.c_str()));
	return list;
}

//---------------------------------------------------------------------------------------------
template<class T>
void median3(vector<T> & ants ,int left ,int right){
	
	int center = (left + right)/2;
	if(ants[center].first < ants[left].first){
		swap(ants[center],ants[left]);
	}
	if(ants[right].first < ants[left].first){
		swap(ants[right],ants[left]);
	}
	if(ants[right].first < ants[center].first){
		swap(ants[center],ants[right]);
	}
	swap(ants[center],ants[right]);
}

template<class T>
int partition(vector<T> &ants,int p ,int r){
	double x = ants[r].first;
	int i = p-1;
	for(int j = p ; j < r ;j++){
		if(ants[j].first <= x){
			i++;
			swap(ants[i],ants[j]);
		}
	}
	swap(ants[i+1],ants[r]);
	return i+1;
}

template<class T>
void quicksort(vector<T> &ants,int p , int r ){
	if(p < r){
		median3(ants,p,r);
		int q = partition(ants,p,r);
		quicksort(ants, p ,q-1);
		quicksort(ants, q+1 ,r);
	}
}

template<class T>
void quicksort(vector<T> &ants){
	quicksort(ants,0,ants.size()-1);
}
//-----------------------------------------------------------------------------------------
template<class T>
void insertionSort(vector<T> &list){

	T key ;
	int id;
	for(unsigned int i=1;i<list.size();i++){
		key = list[i];
		int j = i-1;
		while(j>=0 && list[j].first > key.first){
			list[j+1] = list[j];
			j--;
		}
		list[j+1] = key;
	}
}

int mostSignificant(double n,double decimal){

	int most = int(fmod(n/decimal,10));

	return most;
}
//-------------------------------------------------------------------------
template<class T>
void bucketsort2(vector<T> &ants){
	
	vector<vector<Ant_pair> > buckets(10) ;
	int perBucket;
	double max = ants[0].first;
	double min = ants[0].first;
	for(unsigned int i=1;i<ants.size();i++){
		if(ants[i].first > max){
			max = ants[i].first;
		}
		if(ants[i].first < min){
			min = ants[i].first;
		}
	}
	
	perBucket = int((max-min)/10.0) +1;
	int position;
	for(unsigned int i=0;i<ants.size();i++){
		position = int((ants[i].first-min)/perBucket);
		buckets[position].push_back(ants[i]);
	}
	
	for(unsigned int i=0 ; i<buckets.size();i++){
		insertionSort(buckets[i]);
	}
	
	position = 0;
	for(unsigned int i=0;i<buckets.size();i++){
		for(unsigned int j=0;j<buckets[i].size();j++){
			ants[position++] = buckets[i][j];
		}
	}
	
}
//-----------------------------------------------------------------------------
template<class T>
void bucketsort(vector<T> &ants){
	vector<vector<Ant_pair> > buckets(10) ;
	
	for(unsigned int i=0;i<ants.size();i++){
		if(ants[i].first == 1){
			buckets[9].push_back((ants[i]));
		}
		else{
			buckets[mostSignificant(ants[i].first,pow10(-1))].push_back((ants[i]));
		}
	}
	
	for(unsigned int i=0 ; i<buckets.size();i++){
		insertionSort(buckets[i]);
	}
	
	int position = 0;
	for(unsigned int i=0;i<buckets.size();i++){
		for(unsigned int j=0;j<buckets[i].size();j++){
			ants[position] = buckets[i][j]  ;
			position++;
		}
	}
}


template<class T>
void countsort(vector<T> &ants ,double dec){
	vector<int> 	C(10);
	vector< Ant_pair > B(ants.size());		
	vector<int> 	A(ants.size());
	
	for(unsigned int i=1;i<ants.size();i++){
		A[i] = mostSignificant(pow10(5)*ants[i].first,dec);
	}
	for(unsigned int i=1;i<ants.size();i++){
		C[A[i]] = C[A[i]] +1;
	}
	for(int i=1;i<10;i++){
		C[i] = C[i] + C[i-1];
	}
	for(int i=ants.size()-1;i>=1 ; i--){
		B[C[A[i]]] = ants[i];
		C[A[i]] = C[A[i]] -1; 
	}
	for(unsigned int i=1;i<ants.size();i++){
		ants[i] = B[i];
	}
}

template<class T>
void radixsort(vector<T> &ants){
	
	double max = ants[0].first;
	
	for(unsigned int i=1;i<ants.size();i++){
		if(ants[i].first>max)
			max = ants[i].first;
	}

	int dec = 0;
	while(1){
		if(max < pow10(dec)){
			break;
		}
		dec++;
	}
	for(int i=0;i<dec+5;i++){
		countsort(ants,pow10(i));
	}
}

template<class T>
void shellsort(vector<T> &ants){
	int length = ants.size();
	int h = 1;
	while(h<length/2){
		h = 2*h+1;
	}
	while(h>0){
		for(int i=h;i<length;i++){
			for(int j=i;j>=h && ants[j].first < ants[j-h].first ; j=j-h){
				swap(ants[j],ants[j-h]);
			}
		}
		h = h/2;
	}
}


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

vector<Ant_pair> inputGenerator(int n,int type){
	vector<Ant_pair> input;
	if(type == 0 ){  // same
		for(int i=0;i<n;i++){
			Ant_pair tmp(9987,1);
			input.push_back(tmp);
			
		}
	}
	else if(type == 1 ){  // order
		for(int i=0;i<n;i++){
			Ant_pair tmp(i,1);
			input.push_back(tmp);
			
		}
	}
	else if(type == 2 ){  // random
		for(int i=0;i<n;i++){
			Ant_pair tmp(rand()%10000,1);
			input.push_back(tmp);
			
		}
	}
	
	return input;	
}

int main(int argc,char *argv[]) {
	ifstream inputFile;
	inputFile.open(("hw2.inp"));	
//	inputFile.open(argv[1]);	

	string temp = "";
	int lineNumber = 1;
	int algorithNo ;

	vector<Ant> ants ;
	struct timeval tv_begin;
	struct timeval tv_end;
	int passed_milliseconds;
	
	if(inputFile.is_open()){
		while( !inputFile.eof()){
			getline(inputFile,temp);

			if(lineNumber == 1){
				algorithNo = atoi(temp.c_str());
				lineNumber++;
			}
			else if(lineNumber == 2){
				numberOfAnt = atoi(temp.c_str());
				lineNumber++;
			}
			else{
				vector<double> list = splitString(temp);
				Ant temp ;
				temp.id = lineNumber-2;
				lineNumber++;
				for(int i=0 ; i<N ; i++){
					double t = list[i];
					temp.speed.push_back(t);
				}
				for(int i=N ; i<2*N ; i++){
					double t = list[i];
					temp.time.push_back(t);
				}
				ants.push_back(temp);
				if(temp.id == numberOfAnt)
					break;
			}
		}
		inputFile.close();
	}
	else{
		cout<< "File cannot open" << endl;
	}
	
	//---------------------------------------------------- quick sort  ------------------------------------------------------------
	if(algorithNo == 1 ){	
		ofstream output ;
		vector< Ant_pair > timeAndId;
		output.open("hw2.out");

		for(int j = 0; j<N ; j++){
			timeAndId.clear();
			for(int i=0 ; i < numberOfAnt ; i++ ){
				double speed = ants[i].speed[j];
				double time = 100.0 / speed;
				ants[i].time[j] +=  time;
				Ant_pair p(ants[i].time[j],i+1);
				timeAndId.push_back(p);
			}
			gettimeofday(&tv_begin, NULL);
			quicksort(timeAndId);
			gettimeofday(&tv_end, NULL);
			passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
			output << passed_milliseconds << endl;
			// set points
			double isEqual = -1;
			int k = -1;
			for(unsigned int i= 0 ; i<ants.size();i++){
				if(isEqual != timeAndId[i].first){
					k++;
				}
				isEqual = timeAndId[i].first ;
				ants[timeAndId[i].second-1].point += numberOfAnt-k;
			}
		}
		vector<Ant_pair> pointVector;
		for(int i=0;i<numberOfAnt;i++){
			Ant_pair p(ants[timeAndId[i].second-1].point,timeAndId[i].second);
			pointVector.push_back(p);
		}
		gettimeofday(&tv_begin, NULL);
		quicksort(pointVector);
		
		gettimeofday(&tv_end, NULL);
		passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
		output << passed_milliseconds ;
		
		for(int i=pointVector.size()-1;i>=0;i--){
			output << endl;
			output << pointVector[i].first << " " << pointVector[i].second ;
		}
	}
	else if(algorithNo == 2){	
		ofstream output ;
		vector< Ant_pair > timeAndId;
		output.open("hw2.out");
		
		for(int j = 0; j<N ; j++){
			timeAndId.clear();
			for(int i=0 ; i < numberOfAnt ; i++ ){
				double speed = ants[i].speed[j];
				double time = 100.0 / speed;
				ants[i].time[j] +=  time;
				Ant_pair p(ants[i].time[j],i+1);
				timeAndId.push_back(p);
			}
			gettimeofday(&tv_begin, NULL);
			shellsort(timeAndId);
			gettimeofday(&tv_end, NULL);
			passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
			output << passed_milliseconds << endl;
			// set points
			double isEqual = -1;
			int k = -1;
			for(unsigned int i= 0 ; i<ants.size();i++){
				if(isEqual != timeAndId[i].first){
					k++;
				}
				isEqual = timeAndId[i].first ;
				ants[timeAndId[i].second-1].point += numberOfAnt-k;
			}
		}
		vector< Ant_pair > pointVector;
		for(int i=0;i<numberOfAnt;i++){
			Ant_pair p(ants[timeAndId[i].second-1].point,timeAndId[i].second);
			pointVector.push_back(p);
		}
		gettimeofday(&tv_begin, NULL);
		shellsort(pointVector);
		gettimeofday(&tv_end, NULL);
		passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
		output << passed_milliseconds ;
		
		for(int i=pointVector.size()-1;i>=0;i--){
			output << endl;
			output << pointVector[i].first << " " << pointVector[i].second ;
		}
	}

	//---------------------------------------------------- bucket sort  ------------------------------------------------------------
	else if(algorithNo == 3){	// bucket sort
		vector< Ant_pair > timeAndId;
		ofstream output ;
		output.open("hw2.out");
		for(int j=0;j<N;j++){
			timeAndId.clear();
			for(int i = 0;i<numberOfAnt;i++){
				double time1 = 100.0 / ants[i].speed[j];
				time1 += ants[i].time[j];
				double averageSpeed = 100.0/time1;
				Ant_pair p(averageSpeed,i+1);
				timeAndId.push_back(p);
			}
			
			gettimeofday(&tv_begin, NULL);
			bucketsort(timeAndId);
			gettimeofday(&tv_end, NULL);
			passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
			output << passed_milliseconds << endl;
			
			double isEqual = -1;
			int k=-1;
			for(int i=ants.size()-1; i>=0 ; i--){
				
				if(isEqual != timeAndId[i].first){
					k++;
				}
				isEqual  = timeAndId[i].first;
				ants[timeAndId[i].second-1].point += numberOfAnt-k;
			}
		}
		vector< Ant_pair > pointVector;
		
		for(int i=0;i<numberOfAnt;i++){
			Ant_pair p(ants[timeAndId[i].second-1].point,timeAndId[i].second);
			pointVector.push_back(p);
		}
		
		gettimeofday(&tv_begin, NULL);
//		bucketsort(pointVector);
		bucketsort2(pointVector);
		gettimeofday(&tv_end, NULL);
		passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
		output << passed_milliseconds ;
		for(int i= pointVector.size()-1;i>=0;i--){
			output << endl;
			output << pointVector[i].first << " " << pointVector[i].second ;
		}
		
	}
	//---------------------------------------------------- radix sort -----------------------------------------------------------------------------
	else if(algorithNo == 4){	// radix sort
		ofstream output ;
		output.open("hw2.out");
		
		vector< Ant_pair > timeAndId;
		
		for(int j = 0; j<N ; j++){
			
			timeAndId.clear();
			Ant_pair tmp(-1,-1);
			timeAndId.push_back(tmp);
			
			for(int i=0 ; i < numberOfAnt ; i++ ){
				double speed = ants[i].speed[j];
				double time = 100.0 / speed;
				ants[i].time[j] +=  time;
				Ant_pair tmp(ants[i].time[j],i+1);
				timeAndId.push_back(tmp);
				
			}
			gettimeofday(&tv_begin, NULL);
			radixsort(timeAndId);
			gettimeofday(&tv_end, NULL);
			passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
			output << passed_milliseconds << endl ;
			double isEqual = -1;
			int k = -1;
			for(unsigned int i= 1 ; i<timeAndId.size();i++){
				if(isEqual != timeAndId[i].first){
					k++;
				}
				isEqual = timeAndId[i].first ;
				ants[timeAndId[i].second-1].point += numberOfAnt-k;
			}
		}
		vector< pair<double,int> > pointVector;
		Ant_pair tmp(-1,-1) ;
		pointVector.push_back(tmp);
		for(unsigned int i=1;i<timeAndId.size();i++){
			Ant_pair tmp;
			tmp.first = ants[timeAndId[i].second-1].point;
			tmp.second = timeAndId[i].second;
			pointVector.push_back(tmp);
		}
		gettimeofday(&tv_begin, NULL);
		radixsort(pointVector);
		gettimeofday(&tv_end, NULL);
		passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
		output << passed_milliseconds ;
		for(int i=pointVector.size()-1;i>=1;i--){
			output << endl;
			output << pointVector[i].first << " " << pointVector[i].second;
		}		
	}	
	vector<string> names ;
	names.push_back("same");
	names.push_back("ordered");
	names.push_back("random");
	vector<int> inputsize ;
	inputsize.push_back(50);
	inputsize.push_back(100);
	inputsize.push_back(500);
	inputsize.push_back(1000);
	inputsize.push_back(5000);
	inputsize.push_back(10000);
	
	for(int j=0;j<3;j++){
		cout << names[j] << ":\t";
		if(j!=1){
			cout << "\t" ;
		}
		for(int i=0;i<inputsize.size();i++){
			vector<Ant_pair> input = inputGenerator(inputsize[i],j);
			gettimeofday(&tv_begin, NULL);
			bucketsort2(input);
			gettimeofday(&tv_end, NULL);
			passed_milliseconds = (tv_end.tv_sec - tv_begin.tv_sec)*1000 + (tv_end.tv_usec - tv_begin.tv_usec)/1000;
			cout << passed_milliseconds << "\t" ;
		}	
		cout << endl;
	}

	return 0;
}

