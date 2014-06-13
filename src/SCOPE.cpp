
#include <algorithm>
#include <sys/stat.h>
#include <sys/wait.h>
#include <string>
#include <thread>
#include <mutex>
#include <atomic>
#include "SCOPE.h"
using namespace std;

//#define numTrain 1000

//Global variables
string seq_file = "";                //FASTA/FASTQ input
string out_file= "";          //FASTA output
string file_type = "illumina";   //File type
string file_format = "fasta";
string out_format = "fasta";
string polyType = "A";     //Type of homopolymer
bool doTrim = false;
bool wantDetails = false;
bool printEverything = false;
bool printBestAlign = false;
int filter_width = 12;
int minLength = 10;
int trashLength = 20;
int numTrain = 1000;
int lapK = 1;
int edgeB = 4;
int numStates = 2;
int percentFilter = 65;//The identity filter
bool zeroBased = false;
int right_gap = 30;  //Distance from ends to the beginning of the poly(A)/poly(T) tail
int left_gap = 30;   //Distance from the ends to the end of the poly(A)/poly(T) tail
//  left gap .... AAAAAAA .... right gap
bool identifyTails = false;     //Identify poly(A)/poly(T) tails?
bool polyFront = false;
bool polyBack = false;
bool default_parameters = true; //Determines if default parameters should be used
bool retrain = true;
int numThreads = 1;
int samplingFrequency = 1;
bool mask = false;
vector<string> trainSeqs;
vector<string> backTrainSeqs;//background training seqs
//Files to keep track of stats
//Basic homopolymer detection


std::mutex filelock;

/*This file contains
# poly(A) sequences
# other sequences
# trashed sequences*/
//For basic homopolymer detection
std::atomic<int> polys = ATOMIC_VAR_INIT(0);
std::atomic<int> rejects = ATOMIC_VAR_INIT(0);
std::atomic<int> trashed = ATOMIC_VAR_INIT(0);
//For trimming process
std::atomic<int> polyT_ = ATOMIC_VAR_INIT(0);
std::atomic<int> poly_T = ATOMIC_VAR_INIT(0);
std::atomic<int> polyA_ = ATOMIC_VAR_INIT(0);
std::atomic<int> poly_A = ATOMIC_VAR_INIT(0);
std::atomic<int> polyTT = ATOMIC_VAR_INIT(0);
std::atomic<int> polyAA = ATOMIC_VAR_INIT(0);
std::atomic<int> polyTA = ATOMIC_VAR_INIT(0);
std::atomic<int> polyAT = ATOMIC_VAR_INIT(0);
std::atomic<int> totalSeqs = ATOMIC_VAR_INIT(0);
string stats_file;

#define MAX_LOOP   1000  //Maxinum number of times that window training can scan through a file

void printArgs(){
	cout<<"Input file: "<<seq_file<<endl;                //FASTA/FASTQ input
	cout<<"Output file: "<<out_file<<endl;                //FASTA output
	cout<<"File type: "<<file_type<<endl;   //File type
	cout<<"Input File Format: "<<file_format<<endl;
	cout<<"Output File Format: "<<out_format<<endl;
	cout<<"polyType: "<<polyType<<endl;     //Type of homopolymer
	cout<<"Filter Width: "<<filter_width<<endl;
	cout<<"Edge MinLength: "<<edgeB<<endl;
	cout<<"Boundary States: "<<numStates<<endl;
	cout<<"Mininum Length: "<<minLength<<endl;
	cout<<"Maxinum Training Set: "<<numTrain<<endl;
	cout<<"Laplacian Smoothing Parameter: "<<lapK<<endl;
	cout<<"Details: "<<wantDetails<<endl;
	cout<<"Zero Based: "<<zeroBased<<endl;
	cout<<"Print Everything: "<<printEverything<<endl;
	cout<<"Print Best Alignment: "<<printBestAlign<<endl;

}

static char vDetails[] = "\
SCOPE: Sequence Classification Of homoPolymer Emissions\n\
This software is designed to accurately identify\n\
homopolymers with cDNA sequences using Hidden Markov Models\n\
Miami University, Oxford, Ohio\n\
Bug Report:\n\
Version 0.21: Nov-23-2012";

static char CL_ARGS[] = "\
Command line parameters: \n\
    Input:\n\
       -i [input file] (required) \n\
          the fastq input file or the fasta input file\n\
    Output:\n\
       -o [output file](required)\n\
       --out_format [output file format] \n\
             (default = fasta)\n\
             fasta or fastq\n\
       --details [output details]\n\
             outputs more information including alignment scores,\n\
             homopolymer length, and percent identity\n\
       -z [zero index] \n\
              Output format is printed in zero based indexing, half open intervals \n\
              By default it is printed in one based indexing, closed intervals\n\
    Search Type:\n\
       --homopolymer_type homopolymer type [N|A|G|C|T]\n\
             e.g. option A is a polyA tail\n\
             (default = A)\n\
       --poly\n\
             searches for homopolymers only in the front and back ends\n\
       --trim\n\
             trims poly(A)/poly(T) tails\n\
       --mask\n\
             masks poly(A)/poly(T) tails\n\
    Tool parameters:\n\
       --filter_width filter width\n\
             Size of the sliding window\n\
             (default = 8 base pairs) \n\
       --minLength Mininum homopolymer Length \n\
             (default = 10 base pairs) \n\
       --minIdentity  The minimum identity a homopolymer can have\n\
             (default= 80)\n\
       --edge_minLength Edge boundary MinLength\n\
          (default=6)\n\
       --edge_states Number of states at boundaries\n\
          (default = 1)\n\
       --sampling_frequency  determines how often sequences should be sampled for training\n\
          (default = 1)\n\
       --numTrain Number of training sequences\n\
             (default = 1000) \n\
       --left_gap Distance minLength of beginning of the poly(A)/poly(T) to read end\n\
       --right_gap  Distance minLength of end of the poly(A)/poly(T) to read end\n\
       --no_retrain Disables Baum Welch training\n\
       --numThreads\n\
    Help:\n\
       --help help\n\
       --version version information";
static int verbose_flag;

void parse_cmdline(int argc, char** argv){
	opterr = 0;

	//seq_file;                //FASTA/FASTQ input
	//out_file;                //FASTA output
	file_type = "illumina";   //File format
	polyType = "A";     //Type of homopolymer
	lapK = 1;                //Laplacian smoothing constant
	filter_width = 12;
	numTrain =1000;
	out_format = "fasta";
	edgeB = 4;
	left_gap = 30;
	right_gap = 30;
	int c;
	//numThreads = 1;
	while (1)
	{
		static struct option long_options[] =
		{
				/* These options set a flag. */
				{"verbose", no_argument,       &verbose_flag, 1},
				{"brief",   no_argument,       &verbose_flag, 0},
				/* These options don't set a flag.
                  We distinguish them by their indices. */
				{"i",required_argument,0, 'i'},
				{"file_type",required_argument,0, 's'},
				{"o",required_argument,0, 'o'},
				{"out_format",required_argument,0, 'r'},
				{"print_all",no_argument,0, 'e'},
				{"print_best_align",no_argument,0, 'a'},
				{"details",no_argument,0, 'd'},
				{"z",no_argument,0, 'z'},
				{"homopolymer_type",required_argument,0, 't'},
				{"poly",no_argument,0,'p'},
				{"trim",no_argument,0,'T'},
				{"filter_width",required_argument,0, 'w'},
				{"minLength",required_argument,0, 'm'},
				{"minIdentity",required_argument,0, 'I'},
				{"numTrain",required_argument,0,'n'},
				{"edge_minLength",  required_argument,0, 'b'},
				{"edge_states",  required_argument, 0, 'x'},
				{"left_gap",  required_argument, 0, 'l'},
				{"right_gap",  required_argument, 0, 'g'},
				{"no_retrain", no_argument,0,'V'},
				{"sampling_frequency", required_argument,0,'S'},
				{"numThreads", required_argument,0,'N'},
				{"version",  required_argument, 0, 'v'},
				{"help",    required_argument, 0, 'h'},
				{"masks",  no_argument, 0, 'j'},
				{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "i:o:s:f:r:deat:N:pw:m:I:n:l:g:S:k:b:x:zVjh::v::D::",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{

		case 'i':seq_file = optarg;break;
		case 'o':out_file = optarg;break;
		case 's':file_type = optarg;break;
		case 'r':
			out_format = optarg;
			if(out_format!="fasta" and out_format!="fastq"){
				cout<<"Need a fasta or fastq file"<<endl;
				exit(0);}
			break;
		case 'd':{wantDetails = true;break;}
		case 'e':{printEverything = true;break;}
		case 'a':{printBestAlign = true;break;}
		case 'N':{numThreads = atoi(optarg);break;}
		case 'T':{
			doTrim = true;
			identifyTails = true;
			default_parameters = false;
			break;}
		case 't':{
			string s = optarg;
			if(s=="A"||s=="G"||s=="C"||s=="T"||
					s=="N"||s=="AT"||s=="CG")
				polyType = optarg;
			else{
				cout<<"Unrecognized homopolymer"<<endl;
				exit(1);}
			break;}
		case 'p':identifyTails = true;default_parameters = false;break;
		case 'w':filter_width = atoi(optarg);default_parameters = false;break;
		case 'm':minLength = atoi(optarg);default_parameters = false;break;
		case 'I':percentFilter = atoi(optarg);default_parameters = false;break;
		case 'n':numTrain = atoi(optarg);default_parameters = false;break;
		case 'l':left_gap = atoi(optarg);identifyTails = true;default_parameters = false;break;
		case 'g':right_gap = atoi(optarg);identifyTails = true;default_parameters = false;break;
		case 'k':lapK = atoi(optarg);default_parameters = false;break;
		case 'b':edgeB = atoi(optarg);default_parameters = false;break;
		case 'x':numStates = atoi(optarg);default_parameters = false;break;
		case 'z':zeroBased=true;break;//zeroBased = ((string)optarg=="T");break;
		case 'S':samplingFrequency=atoi(optarg);break;
		case 'V':retrain=false;break;
		case 'j':mask = true; identifyTails = true;break;
		case 'v':
			printf("%s", vDetails);
			cout<<endl;
			//abort();
			exit(-1);
			break;
		case 'h':
			printf("%s", CL_ARGS);
			cout<<endl;
			//abort();
			exit(-1);
			break;
		case 'D':
			seq_file = "test1.fa";
			out_file = "test1_out.fa";
			printBestAlign = true;
			polyType = "AT";
			//file_type = "illumina";
			//file_format = "fastq";
			printEverything = false;
			numStates = 1;
			break;
		case '?':
			if (optopt == 'i' ||
					optopt == 'o')

				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
						"Unknown option character `\\x%x'.\n",
						optopt);
			exit(1);
			break;
		default:
			printf("%s", CL_ARGS);
			abort ();
			break;
		}
	}

	/* Instead of reporting \91--verbose\92
    and \91--brief\92 as they are encountered,
    we report the final status resulting from them. */
	if (verbose_flag)
		puts ("verbose flag is set");

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}

	//       exit (0);

}

void load_defaults(){
	if(file_type=="illumina"){
		filter_width = 12; minLength = 7; lapK = 1;
		edgeB = 1; numStates = 1;right_gap = 25;left_gap = 25;
	}else if(file_type=="sanger"){
		filter_width = 15; minLength = 8;lapK = 1; edgeB = 2;
		numStates = 1; right_gap = 25; left_gap = 25;
	}else{
		minLength = 6; filter_width = 13;  lapK = 1;
		edgeB = 0; numStates = 1; right_gap = 25; left_gap = 25;
	}
}

//std::string readline(ifstream& in){
//  std::string line;
//  std::getline(in,line);
//
//}

bool writeFastaBlock(vector<string> read,ostream& o){
	if(read[0]=="END_OF_FILE" or read[0]=="%.$?/&?&),%")
		return 1;//reached the end
	int c = 60;
	int i = 0;
	string s = read[1];
	while(i<(int)s.length()){
		string str1 = s.substr(0,i);
		string str2 = s.substr(i,s.length());
		s = str1+'\n'+str2;
		i+=c+1;
	}
	if(read[0][0]=='@')
		read[0][0] = '>';
	o<<read[0];
	o<<s<<endl;
	return 0;//Successfully wrote
}

//Reads in one fasta sequence at a time
//pair<string,string> = <description,sequence>
vector<string> readFastaBlock(istream& f){
	vector<string> pair;
	if(f.eof()){
		pair.push_back("END_OF_FILE");
		pair.push_back("END_OF_FILE");
		return pair;
	}else{
		string tmp;
		string cur_des;
		string cur_seq;
		getline(f,cur_des);
		getline(f,tmp);
		//std::cout<<"Description: "<<cur_des<<std::endl;
		cur_seq+=tmp;
		//std::cout<<tmp<<std::endl;
		while(!f.eof() and f.peek()!='>'){
			while(tmp=="\n");
			getline(f,tmp);
			cur_seq+=tmp;
			//std::cout<<tmp<<std::endl;
		}
		

		pair.push_back(cur_des);
		pair.push_back(cur_seq);
		return pair;
	}
}

//quadruple<string,string,string,string> = <description, sequence, description, quality>
vector<string> readFastqBlock(istream& f){
	vector<string> quadruple;
	if(f.eof()){
		return vector<string>(4,"END_OF_FILE");
	}
	string tmp;
	string block;
	getline(f,tmp); quadruple.push_back(tmp);//cout<<tmp<<endl;
	getline(f,tmp); quadruple.push_back(tmp);//cout<<tmp<<endl;
	getline(f,tmp); quadruple.push_back(tmp);//cout<<tmp<<endl;
	getline(f,tmp); quadruple.push_back(tmp);//cout<<tmp<<endl;
	//tmp<<f;quadruple.push_back(tmp);
	//tmp<<f;quadruple.push_back(tmp);
	//tmp<<f;quadruple.push_back(tmp);
	//tmp<<f;quadruple.push_back(tmp);
	//quadruple.push_back(block);
	if(quadruple[0]=="" or
			quadruple[1]=="" or
			quadruple[2]=="" or
			quadruple[3]=="")
		return vector<string>(4,"END_OF_FILE");
	return quadruple;
}

bool writeFastqBlock(vector<string> s,ostream& o){
	if(s[0]=="END_OF_FILE" or s[0]=="%.$?/&?&),%")
		return 1;//reached the end
	o<<s[0]<<endl;
	o<<s[1]<<endl;
	o<<s[2]<<endl;
	o<<s[3]<<endl;
	return 0;
}
/*Get the title immediately following the >*/
std::string getTitle(std::string s){
	std::stringstream ss;
	ss<<s;
	return ss.str();
}

bool isLower(char s){
	if((int)s>=97 and (int)s<=122)
		return true;
	else
		return false;
}
string toUpper(string s){
	for(int i = 0; i<(int)s.length(); i++){
		s[i] = s[i]-32;
	}
	return s;
}
string toLower(string s){
	for(int i = 0; i<(int)s.length(); i++){
		s[i] = s[i]+32;
	}
	return s;
}
string setMask(string s, int start, int end){//Set to lower case
	for(int i = start; i<end; i++){
		if(s[i]>=97 and s[i]<=122){//is lowercase
			continue;
		}else{
			s[i] = s[i]+32;
		}
	}
	return s;
}
char offset(char s){
	char i = s;
	i = s+32;
	return i;
}

bool SCOPE::baumWelchTrain(build_ghmm& model, std::vector<string>& training_seqs){
	Viterbi mle(model.polyA,
			model.backg,
			model.tranP,
			model.start,
			minLength,
			file_type,
			model.poly_filt[0]);
	pair<Matrix<double>,double> F,B;
	Matrix<double> FM,BM;
	double P;
	//model.print_parameters();
	bool doneTraining = false;
	for(unsigned int j = 0 ; j<training_seqs.size();j++){
		F = mle.forward(trainSeqs[j]);
		B = mle.backward(trainSeqs[j]);

		FM = F.first;
		BM = B.first;
		P = F.second;
		doneTraining = model.baumwelchTrain(FM,BM,P,training_seqs[j]);
		//model.print_parameters();

	}
	return doneTraining;

}
build_ghmm SCOPE::buildHMM(char polymerType){
	//cout<<"Building HMM"<<endl;
	build_ghmm model("illumina",//file_type,
  					 filter_width,
					 polymerType,
					 lapK,
					 edgeB,
					 numStates);
	ifstream fin((char*)seq_file.c_str());
	int B = 0;//count of background seqs
	int A = 0;//count of homopolymer seqs
	//int total = 0;//Total number of sequences
	ofstream out;
	srand(time(NULL));
	//int r = 100;//Sample 1 out of 100 randomly
	while(A<numTrain){
		vector<string> read;
		string s;
		if(file_format=="fasta")
			read = readFastaBlock(fin);
		else
			read = readFastqBlock(fin);
		s = read[1];
		if(isLower(s[0]))
			s = toUpper(s);

		if(s!="END_OF_FILE" and s!="%.$?/&?&),%"){
			if(rand()%samplingFrequency!=0){continue;}
			int detected;
			detected=model.train(s);
			// cout<<A<<" "<<B<<" "<<s<<endl;
			// model.print_parameters();

			if(detected==0 and B<numTrain){
				backTrainSeqs.push_back(s);B++;
			}else if(detected==1){
				trainSeqs.push_back(s);A++;
			}
		}else{break;}
	}

	model.finalize_parameters();
	// if(!didDetect){
	//   return build_ghmm();
	// }
	if(retrain and A>0)
		while(!SCOPE::baumWelchTrain(model,trainSeqs));//Train until convergence
	trainSeqs.clear();
	fin.close();
	fin.clear();
	fin.seekg(0);

	//model.print_parameters();
	return model;
}

string SCOPE::trim(vector<Alignment>& alignment, string obs){
	if(alignment.size()==0)
		return obs;
	int start = 0;
	int end = obs.length();
	Alignment first = alignment[0];
	Alignment last = alignment[alignment.size()-1];
	// int minFront = first.start;
	// int minBack = last.seq.length()-last.end;
	if(!first.isNULL() and !last.isNULL()){
		start = first.end;
		end = last.start;
	}else if(first.isNULL() and !last.isNULL()){
		end = last.start;
	}else if(!first.isNULL() and last.isNULL()){
		start = first.end;
	}
	return obs.substr(start,end-start);

}

string SCOPE::applyMask(vector<Alignment> alignments,
		string obs){
	int b;
	int e;
	string s = obs;
	for(int i = 0; i<(int)alignments.size();i++){
		b = alignments[i].start;
		e = alignments[i].end;
		s = setMask(s,b,e);
	}
	return s;
}

vector<Alignment> SCOPE::viterbiCoord(Viterbi& mle,
		string s){

	vector<int> bC;
	vector<int> eC;
	bool isPoly = false;
	//int i = 0;
	string tmp;
	string d;

	//Viterbi algorithm
	string path = mle.viterbi(s);
	//cout<<"Viterbi state sequence"<<endl;
	//cout<<s<<endl;
	//cout<<path<<endl;
	//cout<<endl;

	int k = 0;
	while(k<(int)path.length()){
		if(path[k]=='B'){//Not in polyA
			if(isPoly==true){
				eC.push_back(k);
				isPoly = false;
			}
		}
		else{//In polyA
			if(isPoly ==false){
				bC.push_back(k);
				isPoly = true;
			}
		}
		k++;
	}
	if(bC.size()>eC.size()){
		//if(isPoly==true){
		eC.push_back(k);}
	//isPoly = false;
	vector<Alignment> alignments;
	for(int i = 0; i<(int)bC.size(); i++){
		Alignment a(s,mle.polyType,
				bC[i],eC[i],zeroBased);
		if(a.percentIdentity()>double(percentFilter)/100.0){
			alignments.push_back(a);
		}
	}
	return alignments;
}


pair<Alignment,Alignment> findBorders(vector<Alignment>& align){
	//If there are no alignments, return (NULL,NULL)
	if(align.size()==0){return make_pair(Alignment(),Alignment());}
	//cout<<"Before ";for(int i = 0 ; i<align.size();i++){cout<<align[i].toString1()<<" ";}cout<<endl;

	/*Sort it so that the first and last homopolymers are the borders*/
	sort(align.begin(),align.end(),[](Alignment x, Alignment y){return x.start<y.start;});
	stable_sort(align.begin(),align.end(),[](Alignment x, Alignment y){return x.end<y.end;});
	//cout<<"After ";for(int i = 0 ; i<align.size();i++){cout<<align[i].toString1()<<" ";}cout<<endl;
	vector<Alignment> result;
	Alignment first = align[0];
	Alignment last = align[align.size()-1];
	int minFront = first.start;
	int minBack = last.seq.length()-last.end;
	if(first.equals(last)){
		if(minFront<minBack){
			last = Alignment();
		}else{
			first = Alignment();
		}
		return make_pair(first,last);
	}else{
		//Check to see if they meet the end conditions, otherwise replace with NULL
		if(minFront>left_gap){first = Alignment();}
		if(minBack>right_gap){last = Alignment();}
		return make_pair(first,last);
	}
}

void SCOPE::formatSequence(vector<Alignment>& alignments,
		string s,string q, string& ns, string& nq){
	//string ns,nq;
	if(doTrim){//If trim, it automatically prints everything
		ns = trim(alignments,s);
		if((int)ns.length()<trashLength){ns = "";nq="";return;}//skip over short sequences
		if(file_format=="fastq"){nq = trim(alignments,q);}
		//if(!printEverything and ns==s){continue;}
	}else{
		ns = applyMask(alignments,s);
		if(file_format=="fastq"){nq = q;}
	}
}


/*Assuming poly(A) or poly(T) and that
  it doesn't have both at once*/
Alignment findClosest(vector<Alignment> align){
	cout<<"Found boundary homopolymers"<<endl;
	//cout<<align.size()<<endl;
	if(align.size()==0)
	{
		Alignment a("",'_',0,0);
		return a;
	}
	Alignment bestAlign = align[0];
	int minFront = (align[0].polyType=='T')?(align[0].start):(align[0].seq.length()-align[0].end);

	for(unsigned int i = 1; i<align.size(); i++){
		int curLen = (align[i].polyType=='T')?(align[i].start):(align[i].seq.length()-align[i].end);
		if(curLen<minFront){
			minFront = curLen;
			bestAlign = align[i];
		}
	}
	//int minEnd = minFront+bestAlign.length;
	int minEnd = (bestAlign.polyType=='T')?(bestAlign.seq.length()-bestAlign.end):(bestAlign.start);

	//cout<<minLen<<endl;
	//cout<<minFront<<" "<<minEnd<<" ";
	if(minFront>left_gap or minEnd<right_gap){
		Alignment a("",'_',0,0);
		//cout<<minFront<<" "<<minEnd<<endl;
		return a;
	}else{
		//cout<<endl;
		return bestAlign;
	}
}

string SCOPE::formatDescription(vector<Alignment>& alignments){
	string d = "";
	if(wantDetails==false){
		if((alignments.size())!=0)
			d+=alignments[0].toString1();
		for(int x = 1; x<(int)alignments.size();x++){
			d+=alignments[x].toString1();}
	}else{
		if((alignments.size())!=0)
			d+=alignments[0].toString2();
		for(int x = 1; x<(int)alignments.size();x++){
			d+=alignments[x].toString2();}
	}
	return d;
}

void SCOPE::formatAlignments(vector<Alignment>& alignments){
	if(alignments.size()==0 and !printEverything){return;}
	pair<Alignment,Alignment> p;
	if(identifyTails){
		p = findBorders(alignments);
		alignments.clear();
		alignments.push_back(p.first);
		alignments.push_back(p.second);
	}
}

void SCOPE::printRead(string d, string ns, string nq, ostream& out){
	std::lock_guard<std::mutex> guard(filelock);
	vector<string> newRead;
	if(out_format=="fasta"){
		newRead.resize(2);
		newRead[0] = d; newRead[1] = ns;
		writeFastaBlock(newRead,out);
	}else{
		newRead.resize(4);
		newRead[0] = d; newRead[1] = ns; newRead[3] = nq;
		writeFastqBlock(newRead,out);
	}
}
void SCOPE::printTab(string d, vector<Alignment>& alns, ostream& out, bool isTrash){
	std::lock_guard<std::mutex> guard(filelock);
	std::string title = getTitle(d);
	for(unsigned int i = 0 ; i<alns.size();i++){
		if(alns[i].isNULL())
			continue;
		int start = alns[i].getStart();
		int end = alns[i].getEnd();
		std::string ptype = "";
		if(isTrash)
			ptype = "trash";
		else
			ptype = alns[i].getType();
		out<<title<<"\t"
		   <<start<<"\t"
		   <<end  <<"\t"
		   <<ptype<<endl;
	}
}

//vector<string> readFile(ifstream &f)

bool SCOPE::trash(string ns){
	if((int)ns.size()<trashLength){
		trashed++;
		return true;
	}
	return false;
}
/*Writes to fasta/fastq*/
void SCOPE::writeRecord(vector<Alignment>& alns,
					   string d,
					   string ns,
					   string nq,
					   ostream& fastaout,
					   ostream& tabout,
					   bool isTrash){

	if(alns.size()==0){rejects++;
	}else if(alns[0].isNULL() and alns[1].polyType=='A'){poly_A++;
	}else if(alns[0].isNULL() and alns[1].polyType=='T'){poly_T++;
	}else if(alns[1].isNULL() and alns[0].polyType=='A'){polyA_++;
	}else if(alns[1].isNULL() and alns[0].polyType=='T'){polyT_++;
	}else if(alns[0].polyType=='T' and alns[1].polyType=='T'){polyTT++;
	}else if(alns[0].polyType=='A' and alns[1].polyType=='A'){polyAA++;
	}else if(alns[0].polyType=='T' and alns[1].polyType=='A'){polyTA++;
	}else if(alns[0].polyType=='A' and alns[1].polyType=='T'){polyAT++;
	}else{rejects++;}

	if(not isTrash)
		printRead(d,ns,nq,fastaout);
	printTab(d,alns,tabout,isTrash);
	totalSeqs++;

}


void SCOPE::detect(ifstream& f,
				   ostream& fastaout,
				   ostream& tabout,
				   Viterbi & mle){
	//bool isPoly = false;
	//int i = 0;
	stringstream ss;
	vector<string> read;
	string tmp,s,d,q;
	bool endFlag = 0;
	vector<int> bC,eC;

	/*Read from file*/
	while(!f.eof() and !endFlag){
		{
			/*Critical section: file access*/
			//cout<<"Critical Section 1"<<endl;
			std::lock_guard<std::mutex> guard(filelock);
			if(file_format=="fasta"){
				read = readFastaBlock(f);
				d = read[0]; s = read[1];}
			else{
				read = readFastqBlock(f);
				d = read[0]; s = read[1]; q = read[3];}
			////////////////////////////////
		}

		if(isLower(s[0]))
			s = toUpper(s);
		string ns,nq;
		vector<Alignment> alignments,tmp;

		tmp = viterbiCoord(mle,s);
		for(int j = 0; j<(int)tmp.size();j++){
			alignments.push_back(tmp[j]);}
		formatAlignments(alignments);
		formatSequence(alignments,s,q,ns,nq);
		//d+=formatDescription(alignments);
		vector<string> newRead = read;
		bool isTrash = SCOPE::trash(ns);
		printTab(d,alignments,tabout,isTrash);
		if(not isTrash){
			if(alignments.size()==0){rejects++;}
			else{polys++;}
			printRead(d,ns,nq,fastaout);
		}
		alignments.clear();
	}
}

void SCOPE::trimProcess(ifstream& in,
						ostream& fastaout,
						ostream& tabout,
						Viterbi & modelA,
						Viterbi & modelT){
	stringstream ss;
	vector<string> read;
	string tmp,s,d,q;
	bool endFlag = 0;
	/*Read from file*/
	while(!in.eof() and !endFlag){
		{
			/*Critical section: file access*/
			std::lock_guard<std::mutex> guard(filelock);
			if(file_format=="fasta"){
				read = readFastaBlock(in);
				d = read[0]; s = read[1];}
			else{
				read = readFastqBlock(in);
				d = read[0]; s = read[1]; q = read[3];}
			////////////////////////////////
		}
		if(isLower(s[0])){s = toUpper(s);}
		string ns,nq;
		vector<Alignment> alignments,polyA,polyT;
		if(!modelA.isNULL()){polyA = viterbiCoord(modelA,s);}
		if(!modelT.isNULL()){polyT = viterbiCoord(modelT,s);}
		for(int j = 0; j<(int)polyA.size();j++){alignments.push_back(polyA[j]);}
		for(int j = 0; j<(int)polyT.size();j++){alignments.push_back(polyT[j]);}
		formatAlignments(alignments);
		formatSequence(alignments,s,q,ns,nq);
		//d+=formatDescription(alignments);
		vector<string> newRead = read;

		bool isTrash = SCOPE::trash(ns);
		SCOPE::writeRecord(alignments,
						   d,ns,nq,
						   fastaout,
						   tabout,
						   isTrash);
		alignments.clear();

	}


}
//Retrieve filename from path
std::string getFileName(std::string s){
	unsigned found_dash = s.find_last_of("/");
	// if(found_dash!=std::string::npos){
	//   return s;
	// }else{
	//   int length = s.length()-found_dash;
	std::string file = s.substr(found_dash+1);
	return file;
	//}
}

int main(int argc, char** argv){
	parse_cmdline(argc,argv);
	printArgs();
	// if(default_parameters)
	//   load_defaults();

	if(out_file==""){
	   cout<<"Need output file!!!"<<endl;
	   exit(-1);
	}
	
	// ostream* out =  new ofstream((char*)out_file.c_str());
	if(seq_file==""){
		cout<<"Need input file!!!"<<endl;
		exit(-1);}
	ifstream* in = new ifstream((char*)seq_file.c_str());
	
	if(in->peek()=='>'){file_format="fasta";}
	else{file_format="fastq";}

	if(in->eof())
		cout<<"Empty input file!!!"<<endl;

	std::string fastaOutFile = out_file;
	std::string tabOutFile = fastaOutFile+std::string(".tab");
	std::ofstream fastaout((char*)fastaOutFile.c_str());
	std::ofstream tabout((char*)tabOutFile.c_str());

	if(doTrim or mask){
		vector<Viterbi> Amodels;//Make copies for each thread
		vector<Viterbi> Tmodels;//Make copies for each thread
		build_ghmm modelT = SCOPE::buildHMM('T');
		build_ghmm modelA = SCOPE::buildHMM('A');
		for(int i = 0 ; i<numThreads; i++){
			Amodels.push_back(Viterbi(modelA.polyA,modelA.backg,modelA.tranP,modelA.start,
					minLength,file_type,modelA.poly_filt[0]));
			Tmodels.push_back(Viterbi(modelT.polyA,modelT.backg,modelT.tranP,modelT.start,
					minLength,file_type,modelT.poly_filt[0]));
		}
		vector<std::thread> threads;
		for(int i = 0; i<numThreads; i++){
			threads.push_back(std::thread(SCOPE::trimProcess,
					std::ref(*in),
					std::ref(fastaout),
					std::ref(tabout),
					std::ref(Amodels[i]),
					std::ref(Tmodels[i])));
		}
		for(int i = 0; i<numThreads; i++){threads[i].join();}

		cout<<"[Total Reads]:"<<totalSeqs<<endl;
		cout<<"[TTTT.................]:"<<polyT_<<"\t"<<(polyT_/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[.................TTTT]:"<<poly_T<<"\t"<<(poly_T/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[AAAA.................]:"<<polyA_<<"\t"<<(polyA_/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[.................AAAA]:"<<poly_A<<"\t"<<(poly_A/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[TTTT.............TTTT]:"<<polyTT<<"\t"<<(polyTT/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[AAAA.............AAAA]:"<<polyAA<<"\t"<<(polyAA/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[TTTT.............AAAA]:"<<polyTA<<"\t"<<(polyTA/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[AAAA.............TTTT]:"<<polyAT<<"\t"<<(polyAT/((double)totalSeqs))*100<<"%"<<endl;
		cout<<"[No homopolymers detected]: "<<rejects<<endl;
		cout<<"[Trashed Reads]: "<<trashed<<endl;
	}else{
	  std::cout<<"Building model"<<std::endl;
	  build_ghmm model = SCOPE::buildHMM(polyType[0]);
	  std::cout<<"model finalized"<<std::endl;
	  if(model.isNULL()){
	    cout<<"No homopolymers detected, try again?"<<endl; exit(-1);
	  }
	  vector<Viterbi> vmodels;//Make copies for each thread
	  for(int i = 0 ; i<numThreads; i++){
	    vmodels.push_back(Viterbi(model.polyA,model.backg,model.tranP,model.start,
				      minLength,file_type,model.poly_filt[0]));}
	  vector<std::thread> threads;
	  for(int i = 0; i<numThreads; i++){
	    threads.push_back(thread(SCOPE::detect,std::ref(*in),
				     std::ref(fastaout),std::ref(tabout),std::ref(vmodels[i])));}
	  for(int i = 0; i<numThreads; i++){threads[i].join();}
	  cout<<"Number of sequences with homopolymers "<<polys<<endl;
	  cout<<"Number of sequences without homopolymers "<<rejects<<endl;
	  cout<<"Number of trashed sequences "<<trashed<<endl;
	}

	return 0;
}

