/*
 *  Copyright(C) 2017, Andrej Prochazka, M. Takechi
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nurex/NucID.h"
#include <regex>
#include <string>
#include <cctype>

namespace nurex{

std::regex regexp_nuc("^(\\d{1,3})?([A-z]{1,3})(\\d{1,3})?\\+?$");

std::array<int,2> nucleus_from_symbol(std::string symbol){
    std::array<int,2>res{0,0};
	std::smatch m;
	/*
	int a_size = 0;
	int z_size = 0;
	for(int i=0;i<symbol.size();i++){
		if(!isdigit(symbol[i])){
			a_size = i;
			break;
		}
	};
	for(int i=a_size;i<symbol.size();i++){
		if(isalpha(symbol[i])){
			z_size++;
		}		
	};
	if(z_size == 0 || z_size>3)return res;
	if(a_size == 0 || a_size>3)return res;
	int z = get_element_z(symbol.substr(a_size,z_size));
	if(z<0)return res;
	int a = std::stoi(symbol.substr(0,a_size));
	res[0] = a;
	res[1] = z;
	return res;
	*/
	
    std::regex_match(symbol,m,regexp_nuc);	
	if(m.size()>2){
		int a = std::stoi(m[1]);
		int z = get_element_z(m[2]);
		res[0] = a;
		res[1] = z;
	}
	if(res[1]<=0){
		res[0] = 0;
		res[1] = 0;
		}
    return res;
}

std::map<int,std::string> map_z_symbol{
	{1,"h"},
	{2,"he"},
	{3,"li"},
	{4,"be"},
	{5,"b"},
	{6,"c"},
	{7,"n"},
	{8,"o"},
	{9,"f"},
	{10,"ne"},
	{11,"na"},
	{12,"mg"},
	{13,"al"},
	{14,"si"},
	{15,"p"},
	{16,"s"},
	{17,"cl"},
	{18,"ar"},
	{19,"k"},
	{20,"ca"},
	{21,"sc"},
	{22,"ti"},
	{23,"v"},
	{24,"cr"},
	{25,"mn"},
	{26,"fe"},
	{27,"co"},
	{28,"ni"},
	{29,"cu"},
	{30,"zn"},
	{31,"ga"},
	{32,"ge"},
	{33,"as"},
	{34,"se"},
	{35,"br"},
	{36,"kr"},
	{37,"rb"},
	{38,"sr"},
	{39,"y"},
	{40,"zr"},
	{41,"nb"},
	{42,"mo"},
	{43,"tc"},
	{44,"ru"},
	{45,"rh"},
	{46,"pd"},
	{47,"ag"},
	{48,"cd"},
	{49,"in"},
	{50,"sn"},
	{51,"sb"},
	{52,"te"},
	{53,"i"},
	{54,"xe"},
	{55,"cs"},
	{56,"ba"},
	{57,"la"},
	{58,"ce"},
	{59,"pr"},
	{60,"nd"},
	{61,"pm"},
	{62,"sm"},
	{63,"eu"},
	{64,"gd"},
	{65,"tb"},
	{66,"dy"},
	{67,"ho"},
	{68,"er"},
	{69,"tm"},
	{70,"yb"},
	{71,"lu"},
	{72,"hf"},
	{73,"ta"},
	{74,"w"},
	{75,"re"},
	{76,"os"},
	{77,"ir"},
	{78,"pt"},
	{79,"au"},
	{80,"hg"},
	{81,"tl"},
	{82,"pb"},
	{83,"bi"},
	{84,"po"},
	{85,"at"},
	{86,"rn"},
	{87,"fr"},
	{88,"ra"},
	{89,"ac"},
	{90,"th"},
	{91,"pa"},
	{92,"u"},
	{93,"np"},
	{94,"pu"},
	{95,"am"},
	{96,"cm"},
	{97,"bk"},
	{98,"cf"},
	{99,"es"},
	{100,"fm"},
	{101,"md"},
	{102,"no"},
	{103,"lr"},
	{104,"rf"},
	{105,"db"},
	{106,"sg"},
	{107,"bh"},
	{108,"hs"},
	{109,"mt"},
	{110,"ds"},
	{111,"rg"},
	{112,"cn"},
	{113,"nh"},
	{114,"fl"},
	{115,"mc"},
	{116,"lv"},
	{117,"ts"},
	{118,"og"}};

std::map<std::string,int> map_symbol_z{
	{"h",1},
	{"he",2},
	{"li",3},
	{"be",4},
	{"b",5},
	{"c",6},
	{"n",7},
	{"o",8},
	{"f",9},
	{"ne",10},
	{"na",11},
	{"mg",12},
	{"al",13},
	{"si",14},
	{"p",15},
	{"s",16},
	{"cl",17},
	{"ar",18},
	{"k",19},
	{"ca",20},
	{"sc",21},
	{"ti",22},
	{"v",23},
	{"cr",24},
	{"mn",25},
	{"fe",26},
	{"co",27},
	{"ni",28},
	{"cu",29},
	{"zn",30},
	{"ga",31},
	{"ge",32},
	{"as",33},
	{"se",34},
	{"br",35},
	{"kr",36},
	{"rb",37},
	{"sr",38},
	{"y",39},
	{"zr",40},
	{"nb",41},
	{"mo",42},
	{"tc",43},
	{"ru",44},
	{"rh",45},
	{"pd",46},
	{"ag",47},
	{"cd",48},
	{"in",49},
	{"sn",50},
	{"sb",51},
	{"te",52},
	{"i",53},
	{"xe",54},
	{"cs",55},
	{"ba",56},
	{"la",57},
	{"ce",58},
	{"pr",59},
	{"nd",60},
	{"pm",61},
	{"sm",62},
	{"eu",63},
	{"gd",64},
	{"tb",65},
	{"dy",66},
	{"ho",67},
	{"er",68},
	{"tm",69},
	{"yb",70},
	{"lu",71},
	{"hf",72},
	{"ta",73},
	{"w",74},
	{"re",75},
	{"os",76},
	{"ir",77},
	{"pt",78},
	{"au",79},
	{"hg",80},
	{"tl",81},
	{"pb",82},
	{"bi",83},
	{"po",84},
	{"at",85},
	{"rn",86},
	{"fr",87},
	{"ra",88},
	{"ac",89},
	{"th",90},
	{"pa",91},
	{"u",92},
	{"np",93},
	{"pu",94},
	{"am",95},
	{"cm",96},
	{"bk",97},
	{"cf",98},
	{"es",99},
	{"fm",100},
	{"md",101},
	{"no",102},
	{"lr",103},
	{"rf",104},
	{"db",105},
	{"sg",106},
	{"bh",107},
	{"hs",108},
	{"mt",109},
	{"ds",110},
	{"rg",111},
	{"cn",112},
	{"nh",113},
	{"fl",114},
	{"mc",115},
	{"lv",116},
	{"ts",117},
	{"og",118}};

std::string get_element_symbol(int z){
    std::string res="";
    if(z>0 && z<119){
        res = map_z_symbol[z];
        res[0] = std::toupper(res[0]);
        }
    return res;
    }

int get_element_z(std::string el){
    int res = -1;
    std::transform(el.begin(), el.end(), el.begin(), ::tolower);     
    if(map_symbol_z.count(el)){
        res = map_symbol_z[el];        
        }
    return res;
    }
}