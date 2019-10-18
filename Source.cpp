#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <iostream>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

////////////////////////////////
//// Global variables //////////
////////////////////////////////
vector<vector<string> > counts;
vector<pair<int, int> > missing;
// convergence threshold
float MAX_CHANGE = 0.009;
int MIN_GOOD_VARS = 0;
int MIN_GOOD_VAR_PERCENT = 90;
int prev_good_var = 0;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node {

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
	Graph_Node(string name, int n, vector<string> vals)
	{
		Node_Name = name;

		nvalues = n;
		values = vals;


	}
	string get_name()
	{
		return Node_Name;
	}
	vector<int> get_children()
	{
		return Children;
	}
	vector<string> get_Parents()
	{
		return Parents;
	}
	vector<float> get_CPT()
	{
		return CPT;
	}
	int get_nvalues()
	{
		return nvalues;
	}
	vector<string> get_values()
	{
		return values;
	}
	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT = new_CPT;
	}
	void set_Parents(vector<string> Parent_Nodes)
	{
		Parents.clear();
		Parents = Parent_Nodes;
	}
	// add another node in a graph as a child of this node
	int add_child(int new_child_index)
	{
		for (int i = 0; i < Children.size(); i++)
		{
			if (Children[i] == new_child_index)
				return 0;
		}
		Children.push_back(new_child_index);
		return 1;
	}
	int get_cpt_size() {
		return CPT.size();
	}

};


// The whole network represted as a list of nodes
class network {

	list <Graph_Node> Pres_Graph;

public:
	int addNode(Graph_Node node)
	{
		Pres_Graph.push_back(node);
		return 0;
	}


	int netSize()
	{
		return Pres_Graph.size();
	}
	// get the index of node with a given name
	int get_index(string val_name)
	{
		list<Graph_Node>::iterator listIt;
		int count = 0;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
		{
			if (listIt->get_name().compare(val_name) == 0)
				return count;
			count++;
		}
		return -1;
	}
	// get the node at nth index
	list<Graph_Node>::iterator get_nth_node(int n)
	{
		list<Graph_Node>::iterator listIt;
		int count = 0;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
		{
			if (count == n)
				return listIt;
			count++;
		}
		return listIt;
	}
	//get the iterator of a node with a given name
	list<Graph_Node>::iterator search_node(string val_name)
	{
		list<Graph_Node>::iterator listIt;
		for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
		{
			if (listIt->get_name().compare(val_name) == 0)
				return listIt;
		}

		cout << "node not found\n";
		return listIt;
	}
	string print_nodes_cpt(network Alarm,int var_ind)
	{
			vector<float> node_cpt = Alarm.get_nth_node(var_ind)->get_CPT();
			string ret = "\ttable";
			for (int i = 0; i < node_cpt.size(); i++) {
				float temp = (float)((int)(node_cpt[i] * 10000)) / (float)10000;
				ret += " " + to_string(temp);
			}
			ret += " ;";
			return ret;
	}
	void initialize_min_good_required_vars(int percent)
	{
		MIN_GOOD_VARS = Pres_Graph.size()*percent/100;
	}
	float evaluate()
    {
		float error = 0;

		ifstream myfile("solved_alarm.bif");
		ifstream sirfile("golden_alarm.bif");
		string myline = "", sirline = "";
		while (myline.substr(0, 11).compare("probability") != 0)
		{
			getline(myfile, myline);
			getline(sirfile, sirline);
		}
		// going through each variable cpt table
		int number_of_nodes = Pres_Graph.size();
		for (int k = 0; k < number_of_nodes; k++)
		{
			getline(myfile, myline); // ignoring the probability line
			getline(sirfile, sirline); // ignoring the probability line
			// now we have table line in string variable
			// skiping initial space
			int myindex = 0, sirindex = 0;
			while (sirline[sirindex] != '0') sirindex++;
			while (myline[myindex] != '0') myindex++;
			// now arrived at first number
			// going through all possible numbers
			int my_string_length = myline.size(), sir_string_length = sirline.size(), loop_count = get_nth_node(k)->get_cpt_size();
			for (int i = myindex, j = sirindex, l=0; l<loop_count; i++, j++, l++)
			{
				string my_number = "", sir_number = "";
				// extracting my number
				while (myline[i] != ' ') { my_number += myline[i]; i++; }
				// extracting sir number
				while (sirline[j] != ' ') { sir_number += sirline[j]; j++; }
				// finding the error
				error += abs(stof(my_number) - stof(sir_number));
			}
			getline(myfile, myline); // at {
			getline(sirfile, sirline); // at {
			getline(myfile, myline); // at prob line
			getline(sirfile, sirline); // at prob line
		}
		myfile.close();
		sirfile.close();
		return error;
	}
};

network read_network()
{
	network Alarm;
	string line;
	int find = 0;
	ifstream myfile("alarm.bif");
	string temp;
	string name;
	vector<string> values;

	if (myfile.is_open())
	{
		while (!myfile.eof())
		{
			stringstream ss;
			getline(myfile, line);


			ss.str(line);
			ss >> temp;


			if (temp.compare("variable") == 0)
			{

				ss >> name;
				getline(myfile, line);

				stringstream ss2;
				ss2.str(line);
				for (int i = 0; i < 4; i++)
				{

					ss2 >> temp;


				}
				values.clear();
				while (temp.compare("};") != 0)
				{
					temp = temp.substr(1, temp.size() - 2);
					values.push_back(temp);

					ss2 >> temp;
				}
				Graph_Node new_node(name, values.size(), values);
				int pos = Alarm.addNode(new_node);


			}
			else if (temp.compare("probability") == 0)
			{

				ss >> temp;
				ss >> temp;

				list<Graph_Node>::iterator listIt;
				list<Graph_Node>::iterator listIt1;
				listIt = Alarm.search_node(temp);
				int index = Alarm.get_index(temp);
				ss >> temp;
				values.clear();
				while (temp.compare(")") != 0)
				{
					listIt1 = Alarm.search_node(temp);
					listIt1->add_child(index);
					values.push_back(temp);

					ss >> temp;

				}
				listIt->set_Parents(values);
				getline(myfile, line);
				stringstream ss2;

				ss2.str(line);
				ss2 >> temp;

				ss2 >> temp;

				vector<float> curr_CPT;
				string::size_type sz;
				while (temp.compare(";") != 0)
				{
					curr_CPT.push_back(atof(temp.c_str()));
					ss2 >> temp;
				}
				listIt->set_CPT(curr_CPT);
			}
			else
			{
			}
		}

		if (find == 1)
			myfile.close();
	}

	return Alarm;
}

/////////////////////////////
//// global variable ////////
/////////////////////////////
network Alarm;
////////////////////////////

int get_ind_value_match(string value, vector<string> values) {
	int ind = 0;
	for (ind = 0; ind < values.size(); ind++) {
		if (value.compare(values[ind]) == 0)
			return ind;
	}
	return -1;
}

int get_CPT_index(vector<int> match_ind, vector<int> indexes) {
	int ind = 0;
	int mult = 1;
	for (int i = indexes.size()-1; i >=0; i--) {
		ind += mult * match_ind[i];
		mult *= Alarm.get_nth_node(indexes[i])->get_nvalues();
	}
 	return ind;
}

vector<float> get_CPT(list<Graph_Node>::iterator node) {
	vector<float> temp;
	vector<string> vals = node->get_values();
	vector<string> pars = node->get_Parents();
	vector<int> indexes;
	indexes.push_back(Alarm.get_index(node->get_name()));
	int num = vals.size();
	for (int i = 0; i < pars.size(); i++) {
		num *= Alarm.search_node(pars[i])->get_nvalues();
		indexes.push_back(Alarm.get_index(pars[i]));
	}
	for (int i = 0; i < num; i++) {
		temp.push_back(1);
	}
	// P(A|B,C) -> temp contains first all values of A for a given B and C. 
	//Then the next value of C with B unchnaged, all A values for this new pair of evidence.
	//Then when all values of C exhausted, then next value of B.
	for (int i = 0; i < counts.size(); i++) {
		vector<int> match_index;
		for (int j = 0; j < indexes.size(); j++) 
			match_index.push_back(get_ind_value_match(counts[i][indexes[j]], Alarm.get_nth_node(indexes[j])->get_values())); 
		int cpt_ind = get_CPT_index(match_index, indexes);
		temp[cpt_ind]++;
	}
	int sum = 0;
	num /= vals.size();
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < vals.size(); j++)
			sum += temp[i+(num*j)];
		for (int j = 0; j < vals.size(); j++)
			temp[i+(num*j)] /= sum;
		sum = 0;
	}
	return temp;
}

// apply kl divergence check
bool kl_divergence(vector<float> old_cpt, vector<float> &new_cpt)				//// ####
{
	int cpt_size = old_cpt.size();
	if(cpt_size != new_cpt.size())	{ cout << "Error: sizes of new calculated cpt does not match with old cpt size.\n"; return false;}
	float divergence_value = 0;
	// traversing through all cpt entry
	for(int i=0; i<cpt_size; i++)
	{
		//cout << old_cpt[i] << " " << new_cpt[i] << " " << log(old_cpt[i]) << " " << log(new_cpt[i]) << endl;
		divergence_value += new_cpt[i]*(log(old_cpt[i]) - log(new_cpt[i]));
	}
	// check if converged
	//cout << divergence_value << endl;
	return (abs(divergence_value) < MAX_CHANGE);
}

bool compute_prob() {										//// ####
	bool is_converged = false;
	int good_var = 0;
	for (int i = 0; i < Alarm.netSize(); i++) {
		vector<float> new_cpt = get_CPT(Alarm.get_nth_node(i));
		// converge only if all variable cpts converge
		if(kl_divergence(Alarm.get_nth_node(i)->get_CPT(), new_cpt) == true)	good_var++;
		Alarm.get_nth_node(i)->set_CPT(new_cpt);
	}
	cout << MIN_GOOD_VARS << endl;
	cout << "good variables: " << good_var << endl;
	cout << endl;
	if(good_var < MIN_GOOD_VARS)	is_converged = false;
	if(good_var >= MIN_GOOD_VARS && good_var < prev_good_var) is_converged = true;
	if(good_var > MIN_GOOD_VARS)
		prev_good_var = good_var;
	return is_converged;
}


// initializing cpts
void initialize_cpts()
{
	// modifing cpt of each node in network
	int network_size = Alarm.netSize();
	for (int i = 0; i < network_size; i++)
	{
		// creating initial cpt with value given in arg
		vector<float> initial_cpt;
		list<Graph_Node>::iterator node_it = Alarm.get_nth_node(i);
		int cpt_size = node_it->get_cpt_size();
		int num_vals = node_it->get_nvalues();
		cpt_size /= num_vals;
		//float initial_prob_val = 1.0 / (node_it->get_nvalues());
		//for (int i = 0; i < cpt_size; i++)	initial_cpt.push_back(initial_prob_val);
		for(int k=0; k<num_vals-1; k++){
			for(int j=0; j<cpt_size;j++){
				float random = (float)(rand()/((double)RAND_MAX * (num_vals-1)));
				initial_cpt.push_back(random);
			}	
		}
		for(int j=0; j<cpt_size;j++){
			float temp = 0;
			for(int i=0; i<num_vals-1;i++){
				temp += initial_cpt[j+(cpt_size * i)];
			}
			initial_cpt.push_back(1-temp);
		}
		Alarm.get_nth_node(i)->set_CPT(initial_cpt);
	}
}


float get_CPT_of_child(int child_ind, int row_num, int ind, int var_num) {
	vector<string> pars = Alarm.get_nth_node(child_ind)->get_Parents();
	vector<int> indexes;
	indexes.push_back(child_ind);
	for (int i = 0; i < pars.size(); i++) {
		int temp = Alarm.get_index(pars[i]);
		indexes.push_back(temp);
	}
	vector<int> match_index;
	for (int i = 0; i < indexes.size(); i++) {
		if (indexes[i] == var_num)
			match_index.push_back(ind);
		else
		match_index.push_back(get_ind_value_match(counts[row_num][indexes[i]], Alarm.get_nth_node(indexes[i])->get_values()));
	}
	int cpt_ind = get_CPT_index(match_index, indexes);
	return (Alarm.get_nth_node(child_ind)->get_CPT())[cpt_ind];
}

// missing vector contains pair of integers where first is the record number and the second is the variable index which is missing
void sample_missing_points() {
	for (int i = 0; i < missing.size(); i++) {
		int r_num = missing[i].first;
		int v_num = missing[i].second;
		vector<string> vals;
		vals = Alarm.get_nth_node(v_num)->get_values();
		int ran = rand()%100;
		if(ran >0){
			vector<int> child;
			child = Alarm.get_nth_node(v_num)->get_children();
			vector<float> sample_prob;
			for (int k = 0; k < vals.size(); k++) {
				float prob = 1;
				for (int j = 0; j < child.size(); j++) 
					prob *= get_CPT_of_child(child[j], r_num, k, v_num);
				prob *= get_CPT_of_child(v_num, r_num, k, v_num);
				sample_prob.push_back(prob);
			}
			float sum = 0;
			for (int k = 0; k < sample_prob.size(); k++)
				sum += sample_prob[k];
			for (int k = 0; k < sample_prob.size(); k++) {
				sample_prob[k] /= sum;
				if (k > 0)
					sample_prob[k] += sample_prob[k - 1];
			}
			float rand_num = (float)(rand() / (double)RAND_MAX);
			for (int k = 0; k < sample_prob.size(); k++) {
				if (rand_num <= sample_prob[k]) {
					counts[r_num][v_num] = vals[k];
					break;
				}
			}
		}
		else{
			int random = rand()%(vals.size());
			counts[r_num][v_num]= vals[random];
		}
	}
}


void parse_data_file(string filename, int num_vars) {
	ifstream infile("records.dat");
	string val;
	int ind = 0;
	int record_num = 0;
	vector<string> temp;
	while (infile >> val) {
		if (val.compare("") == 0)
			break;
		val = val.substr(1, val.size() - 2);
		if (ind == 0)
			temp.clear();
		if (val.compare("?") == 0)
			missing.push_back(make_pair(record_num, ind));
		temp.push_back(val);
		if (ind == num_vars - 1) {
			record_num++;
			//cout << record_num << endl;
			counts.push_back(temp);
		}
		ind = (ind + 1) % num_vars;
	}
	initialize_cpts();
	sample_missing_points();
}

void write_final_output() {
	ifstream infile("alarm.bif");
	ofstream myfile;
	myfile.open("solved_alarm.bif");
	string line = "";
	while (getline(infile, line)) {
		if (line.substr(0, 11).compare("probability") == 0)
			break;
		else
			myfile << line << "\n";
	}
	for (int i = 0; i < Alarm.netSize(); i++) {
		myfile << line << "\n";
		myfile << Alarm.print_nodes_cpt(Alarm,i) << "\n";
		getline(infile, line);
		getline(infile, line);
		myfile << line << '\n';
		getline(infile, line);
	}
	infile.close();
	myfile.close();
}

int main(int argc, char *argv[])
{
	Alarm = read_network();
	Alarm.initialize_min_good_required_vars(MIN_GOOD_VAR_PERCENT);
	string filename = "";
	int num_vars = Alarm.netSize();
	srand(time(0));
	//	create_records_storage();
	parse_data_file(filename, num_vars);
	//appling em algorithm loop
	bool is_converged = false;
	int deleteThis=0;
	while(not is_converged and (deleteThis <= 800))
	{
		is_converged = compute_prob();
		sample_missing_points();
		deleteThis++;
		cout << "iteration no: " << deleteThis << endl;
	}
	write_final_output();
	cout << Alarm.evaluate() << endl;
}
