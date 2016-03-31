#include <vector>
#include <iostream>
#include <math.h>
#include <unordered_map>

using namespace std;
using f_pair = pair<float, float>;
using f_vec = vector<f_pair>;
f_vec points   {
	make_pair(4.0, 5.0),	make_pair(3.0, 8.0),	make_pair(0.0, 3.0),
	make_pair(5.0, 5.0), 	make_pair(3.0, 4.0), 	make_pair(1.0, 9.0),
	make_pair(3.0, 2.0), 	make_pair(4.0, 6.0), 	make_pair(2.0, 8.0)};

template<typename val> 
void print(const f_vec& that) {
	for(const auto pr: that)	cout << pr.first << " " << pr.second << endl;
}

template<typename val>
void print(const vector<val>& that) {
	for(const auto v: that) cout << v << " ";
	cout << endl;
}

template<typename val>
void print(const f_pair that) {
	cout << "{" << that.first << " " << that.second << "} " << endl;
}

float euc_distance(const f_pair p1, const f_pair p2) {
	float var1 = p1.first - p2.first;
	float var2 = p1.second - p2.second;
	return sqrt(var1*var1 + var2*var2);
}

void knn(const int k, vector<float>& k_dis, const f_vec& points) {
	k_dis.clear();
	for(const auto p1: points) {
		vector<float> distance;
		for(const auto p2: points) {
			float d_p1p2 = euc_distance(p1, p2);
			if(d_p1p2 == 0.0)	continue;
			else	distance.push_back(d_p1p2);
		}
		sort(distance.begin(), distance.end());
		k_dis.push_back(distance[k - 1]);
	}
}

float find_knn(const vector<float>& k_dis, 
			const f_vec& points,
			const f_pair pt){
	size_t index = 0;
	for(auto itor = points.cbegin(); itor != points.cend(); ++itor) {
		if(*itor == pt)	return k_dis[index];
		index++;
	}
	return 0.0;
}

float reach_dis(const float k_d, const pair<float, float> p1,
			const pair<float, float> p2) {
	float dis = euc_distance(p1, p2);
	return k_d > dis ? k_d : dis;
}

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

void find_nk(f_vec& nk, const f_pair point,
			const f_vec& points) {
	nk.clear();
	unordered_map<pair<float, float>, float, pairhash> umap;
	vector<float> dis_v;
	for(const auto p1: points) {
		float d_p1 = euc_distance(p1, point);
		if(d_p1 == 0.0)	continue;
		//umap.push_back(d_p1p2, point);
		umap.insert(make_pair(p1, d_p1));
		dis_v.push_back(d_p1);
	}
	sort(dis_v.begin(), dis_v.end());
	for(size_t i = 0; i < 3; ++i){
		float f = dis_v[i];
		for(auto itor = umap.begin(); itor != umap.end(); ++itor) {
			if(itor -> second == f)	{
				nk.push_back(itor -> first);
				umap.erase(itor);
				continue;
			}
		}
	}
	if(dis_v[2] == dis_v[3]) {
		for(auto itor = umap.begin(); itor != umap.end(); ++itor) {
			if(itor -> second == dis_v[2])	nk.push_back(itor -> first);
		}
	}
}

float lrd(const f_vec& nk, const f_pair point, 
			const vector<float>& k_dis) {
	const int len = nk.size();
	float sum = 0.0;
	for(const auto pt: nk) {
		float rd = find_knn(k_dis, points, pt);
		float reach = reach_dis(rd, pt, point);
		sum += reach;
	}
	return len/sum;
}


int main() {
	vector<float> k_dis;
	f_vec nk;
	float sum = 0.0;
	cout << "Answer For Part A" << endl;
	auto v = points;
	knn(3, k_dis, points);
	for(size_t i = 0; i < 9; ++i) {
		cout << "\tFor Point {" << v[i].first << " " << v[i].second << "},";
		cout << "\t" 	<< "3-distance is " << "\t" << "sqrt(" 
						<< k_dis[i]*k_dis[i]  << ")" << endl;;
	}
	cout << endl;

	cout << "Answer For Part B" << endl;
	cout << "\tStep 1: Calculate All lrd" << endl;
	for(auto const point: points) {
		vector<pair<float, float>> new_nk;
		find_nk(new_nk, point, points);
		float lrd_x = lrd(new_nk, point, k_dis);
		cout << "\tFor Point {" << point.first << " " << point.second <<"},";
		cout << "\tlrd = " << "\t" << lrd_x << endl;
	}

	cout << endl << "\tStep 2: Calcualte L.O.F" << endl;

	for(const auto point: points) {
		//print(point);

		nk.clear();
		find_nk(nk, point, points);

		const float len = float(nk.size());
		float sum = 0.0;
		float lrd_A = lrd(nk, point, k_dis);

		//print (nk);
		for(const auto p_B: nk) {
			vector<pair<float, float>> new_nk;
			find_nk(new_nk, p_B, points);
			//cout << lrd(new_nk, p_B, k_dis) << endl;
			sum += lrd(new_nk, p_B, k_dis);
		}
		float result = sum / (len * lrd_A);

		cout << "\tFor Point {" << point.first << " " << point.second <<"},";
		cout << "\tL.O.F = " << result << endl;
		//print(point);
		//cout << "LOF = " << sum / (len * lrd_A) << endl << endl;

	}

	return 0;
}

