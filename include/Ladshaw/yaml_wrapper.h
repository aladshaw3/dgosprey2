/*!
 *  \file yaml_wrapper.h yaml_wrapper.cpp
 *	\brief C++ Wrapper for the C-YAML Library
 *	\details This file holds objects, structures, and functions associatied with using the C-YAML library.
 *			A C++ wrapper has been created for the Kirill Simonov (2006) LibYAML library to more easily
 *			store and query information in yaml style input files. The wrapper uses the C-YAML parser to
 *			identify the file structure and store the read in information from that document into an
 *			object using C++ maps. Those maps are hold information in a series of Key-Value pairs as well
 *			as lists of Key-Value pairs. This allows the user to create well organized input files to
 *			change the behavior of simulations. 
 *
 *			The yaml_wrapper is restricted to the same limitations in the C-YAML source code in terms of
 *			how the documents are allowed to be structured for TOKEN based parsing. C-YAML only recognizes
 *			specific tokens and will only allow a certain level of Sub-Header mapping. Therefore, this
 *			wrapper has the same limitations. Below is an example of acceptable formatting for a C-YAML
 *			document.
 *
 			\#Test input file for YAML and SHARK \n
 
			TestDoc1: &hat \n
			--- \n
			- scenario: \n
			  numvar: 25 \n
			  act_fun: DAVIES \n
			  steadystate: FALSE \n
			  t_out: 1 \n
			  pH: 0 \n
 
			- testblock: \n
			  another: block \n
			  - subblock: \n
			    sub: block \n
			... \n
 
			TestDoc2: *hat \n
			--- \n
			- masterspecies: \n
			  "Cl - (aq)": 0 \n
			  "Na + (aq)": 1 \n
			  "H2O (l)": 2 \n
			  3: NaCl (aq) \n
			... \n
 
			TestDoc3: \n
			--- \n
			apple: red \n
			pear: green \n
 
			- array: \#Block \n
			  banana: yellow \n
			  \#List 1 in array \n
			  - list1: \&a \#also a block \n
			    a: 1 \#key : value \n
			    b: 2 \n
			    c: 3 \n
			\#List 2 in array
			  - list2: *a \n
			    a: 4 \n
			    b: 5 \n
			    c: 6 \n
			... \n
 
			TestDoc4: \n
			--- \n
			- anchor: \&anchor \n
			  stuff: to do \n
			- alias: *anchor \n
			  add: to stuff \n
 
			- list: \n
			  - anchored: \&list_anchor \n
			    info: blah \n
			    atta: boy \n
			  - aliased: *list_anchor \n
			    info: bruh \n
			    atta: ber \n
			... \n
 
			\#WARNING: MAKE SURE FILE DOES NOT CONTAIN TABS!!!
 
			TestDoc5: \n
			--- \n
 
			- grab: *anchor \n
			  add2: more adds \n
			  - listcopy: *a \n
 
			- block: {1: 2, 3: 4} \n
			  still: in block \n
			... \n
 
 *  \note You can view the actual yaml example file in the input_files/SHARK/test_input.yml sub-directory of the project folder.
 *  \author Austin Ladshaw
 *	\date 07/29/2015
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for PhD research in the area
 *             of adsorption and surface science. Copyright (c) 2015, all
 *             rights reserved. This copyright only applies to yaml_wrapper.h
 *			   and yaml_wrapper.cpp. \n
 *				
 *				DISCLAIMER: 
				----------
				Niether Austin Ladshaw, nor the Georgia Institute of Technology, is the author or owner of any C-YAML
				Library or source code. Only the files labeld "yaml_wrapper" were created by Austin Ladshaw. Therefore, any C-YAML
				files distributed will be given a portions copyright under the MIT License (see below). For more information on YAML, 
				go to pyyaml.org/wiki/LibYAML.\n
 
 				Portions copyright 2006 Kirill Simonov \n
 
				The MIT License (MIT) \n
 
				Permission is hereby granted, free of charge, to any person obtaining a copy
				of this software and associated documentation files (the "Software"), to deal
				in the Software without restriction, including without limitation the rights
				to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
				copies of the Software, and to permit persons to whom the Software is
				furnished to do so, subject to the following conditions: \n
 
				The above copyright notice and this permission notice shall be included in
				all copies or substantial portions of the Software. \n
 
				THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
				IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
				FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
				AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
				LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
				OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
				THE SOFTWARE. \n
 */


#include "yaml.h"
#include "error.h"
#include <map>
#include <string>
#include <iostream>
#include <utility>
#include <stdexcept>

#ifndef YAML_WRAPPER_HPP_
#define YAML_WRAPPER_HPP_

typedef enum data_type
{ STRING, BOOLEAN, DOUBLE, INT, UNKNOWN } data_type;			///< Enum for valid data types in ValueTypePair

typedef enum header_state { ANCHOR, ALIAS, NONE } header_state;	///< Enum for state of the headers in the yaml_wrapper

/// Value-Type Pair object to recognize data type of a string that was read
/** C++ Object that creates a pair between a read in value as a string and an enum denoting
	what the data type of that string is. This object is primarily used in the other yaml_wrapper
	objects, but can also be used for any string that you want to parse to identify it's type.
	The supported types are denoted in the data_type enum and can be determined automatically 
	by the findType() function or can be specified by the assertType() function. */
class ValueTypePair
{
public:
	ValueTypePair();										///< Default constructor
	~ValueTypePair();										///< Default destructor
	ValueTypePair(const std::pair<std::string,int> &vt);	///< Constructor by pair
	ValueTypePair(std::string value, int type);				///< Construction by string and int
	ValueTypePair(const ValueTypePair &vt);					///< Copy constructor
	
	ValueTypePair& operator=(const ValueTypePair &vt);		///< Equals operator overload
	
	void editValue(std::string value);					///< Edits value to pair with UNKOWN type
	void editPair(std::string value, int type);			///< Creates a paired Value-Type from the given args
	void findType();									///< Determines the data type of the object
	void assertType(int type);							///< Forces a specific data type
	void DisplayPair();									///< Display the pair information
	
	std::string getString();							///< Returns the value of the pair as a string
	bool getBool();										///< Returns the value of the pair as a bool
	double getDouble();									///< Returns the value of the pair as a double
	int getInt();										///< Returns the value of the pair as an int
	std::string getValue();								///< Returns the value of the pair as it was given
	int getType();										///< Returns the type of the pair
	std::pair<std::string,int> &getPair();				///< Returns reference to the actual pair object
	
protected:
	
private:
	std::pair<std::string,int> Value_Type;				///< pair object holding the Value and Type info
	int type;											///< Type of the value
	
};

/// Key-Value-Type Map object creating a map of the KeyValuePair objects
/** C++ Object that creates a map of the KeyValuePair objects. Functions defined here allow the user to
	iterate through this map, access specific keys in the map, edit values associated with those keys,
	find the data types for the values in those keys, ect. The keys are used as an access operator for
	their corresponding value. As such, each key in the map is required to be unique, but the values are
	allowed to be duplicated. */
class KeyValueMap
{
public:
	KeyValueMap();										///< Default constructor
	~KeyValueMap();										///< Default destructor
	KeyValueMap(const std::map<std::string,std::string> &map);///< Construct from a map of strings
	KeyValueMap(std::string key, std::string value);	///< Construct one element in the map
	KeyValueMap(const KeyValueMap &map);				///< Copy constructor
	
	KeyValueMap& operator=(const KeyValueMap &map);			///< Equals overload
	ValueTypePair& operator[](const std::string key);		///< Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	///< Return the ValueType at the give key
	
	std::map<std::string, ValueTypePair > & getMap();					///< Return a reference to the Key_Value map object
	std::map<std::string, ValueTypePair>::const_iterator end() const;	///< Returns a const iterator pointing to the end of the list
	std::map<std::string, ValueTypePair>::iterator end();				///< Returns an iterator pointing to the end of the list
	std::map<std::string, ValueTypePair>::const_iterator begin() const;	///< Returns a const iterator pointing to the beginning of the list
	std::map<std::string, ValueTypePair>::iterator begin();				///< Returns an iterator pointing to the beginning of the list
	
	void clear();														///< Clears the map
	void addKey(std::string key);										///< Adds a key to the object with a default value
	void editValue4Key(std::string val, std::string key);				///< Edits a given value for a pre-existing key
	void editValue4Key(std::string val, int type, std::string key);		///< Edits a value for a pre-existing key and asserts type
	
	void addPair(std::string key, ValueTypePair val);	///< Adds a pair object to the map
	void addPair(std::string key, std::string val);		///< Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int type); ///< Adds a pair object and asserts a type
	void findType(std::string key);						///< Find what data type the value at the key is
	void assertType(std::string key, int type);			///< Assert the given type at the given key
	void findAllTypes();								///< Find all types for all data in map
	void DisplayMap();									///< Print out the map to console
	
	int size();											///< Returns the size of the map
	std::string getString(std::string key);				///< Retrieve the string at the key
	bool getBool(std::string key);						///< Retrieve the boolean at the key
	double getDouble(std::string key);					///< Retrieve the double at the key
	int getInt(std::string key);						///< Retrieve the int at the key
	std::string getValue(std::string key);				///< Retrieve the value at the key
	int getType(std::string key);						///< Retrieve the type at the key
	ValueTypePair& getPair(std::string key);			///< Retrieve the pair at the key
	
protected:
	
private:
	std::map<std::string, ValueTypePair > Key_Value;		///< Map of Keys and Values paired with types
};

/// Object for the Lowest level of Header for the yaml_wrapper
/** C++ Object for sub-headers in a yaml document. This object contains a KeyValueMap that holds a set of
	key-value pairs for data listed under a sub-header in yaml files. It is the lowest allowable recursion
	of headers in a yaml document and so is the base class for Header and Document, which themselves can
	contain KeyValueMaps as well as maps for other header-like objects. 
 
	SubHeaders are recognized by a unique name and/or alias while being put together in other higher 
	document structures. Additionally, each header or sub-header will have a state to denote whether
	the object is a yaml alias, anchor, or niether. This is used in the yaml documents to ensure that
	aliases for anchors have the correct data moved over into the new structures. */
class SubHeader
{
public:
	SubHeader();										///< Default Constructor
	~SubHeader();										///< Default Destructor
	SubHeader(const SubHeader &subheader);				///< Copy constructor
	SubHeader(const KeyValueMap &map);					///< Construction by existing map
	SubHeader(std::string name);						///< Construction by name only
	SubHeader(std::string name, const KeyValueMap &map);///< Construction by name and map
	
	SubHeader& operator=(const SubHeader &sub);				///< Equals overload
	ValueTypePair& operator[](const std::string key);		///< Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	///< Return the ValueType at the give key
	
	KeyValueMap& getMap();									///< Returns reference to the KeyValueMap object
	
	void clear();																///< Empty out data contents
	void addPair(std::string key, std::string val);								///< Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int type);					///< Adds a pair object and asserts a type
	void setName(std::string name);												///< Sets the name of the subheader
	void setAlias(std::string alias);											///< Set the alias without type specification
	void setAlias(std::string alias, int state);								///< Sets the alias and state of the subheader
	void setNameAliasPair(std::string name, std::string alias, int state);		///< Sets the name and alias of the subheader
	void setState(int state);													///< Sets the state of the subheader
	void DisplayContents();														///< Display the contents of the subheader
	
	std::string getName();								///< Return the name of the subheader
	std::string getAlias();								///< Return the alias of the subheader, if one exists
	bool isAlias();										///< Returns true if subheader is an alias
	bool isAnchor();									///< Returns true if subheader is an anchor
	int getState();										///< Returns the state of the subheader
	
protected:
	KeyValueMap Data_Map;								///< A Map of Keys and Values
	std::string name;									///< Name of the subheader
	std::string alias;									///< Name of the alias for the subheader
	int state;											///< State of the header
	
private:
	
};

/// Object for headers in a yaml document (inherits from SubHeader)
/** C++ Object for headers in a yaml document that is built from the SubHeader object already created.
	The chain of inheritance works in this direction because a Header can have both a map of SubHeaders and
	a map of KeyValuePairs. Therefore, the SubHeader object is actually the more generic form of a header.
 
	Since this object inherits from SubHeader, it has access to all it's protected members, including the
	alias, state, name, and KeyValueMap. Operator overloads and other functions are provided to allow the
	user to query both the KeyValueMap and SubHeader for specific information. The names of the SubHeaders
	are also used as it's keys. Make sure all SubHeader keys are unique to this header. */
class Header : SubHeader
{
public:
	Header();											///< Default Constructor
	~Header();											///< Default Destructor
	Header(const Header &head);							///< Copy constructor
	Header(std::string name);							///< Constructor by header name
	Header(const KeyValueMap &map);						///< Constructor by existing map
	Header(std::string name, const KeyValueMap &map);	///< Constructor by name and map
	Header(std::string key, const SubHeader &sub);		///< Constructor by single subheader object
	
	Header& operator=(const Header &head);					///< Equals overload
	ValueTypePair& operator[](const std::string key);		///< Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	///< Return the ValueType at the given key
	SubHeader& operator()(const std::string key);		    ///< Return the SubHeader reference at the given key
	SubHeader operator()(const std::string key) const;		///< Return the SubHeader at the given key
	
	std::map<std::string, SubHeader> & getSubMap();			///< Return the reference to the SubHeader Map
	KeyValueMap & getDataMap();								///< Return the reference to the KeyValueMap
	SubHeader & getSubHeader(std::string key);				///< Return the subheader at the given key
	std::map<std::string, SubHeader>::const_iterator end() const;	///< Returns a const iterator pointing to the end of the list
	std::map<std::string, SubHeader>::iterator end();				///< Returns an iterator pointing to the end of the list
	std::map<std::string, SubHeader>::const_iterator begin() const;	///< Returns a const iterator pointing to the begining of the list
	std::map<std::string, SubHeader>::iterator begin();				///< Returns an iterator pointing to the begining of the list
	
	void clear();											///< Clear out the SubMap, KeyValueMap, and other info
	void resetKeys();										///< Reset the keys of the SubMap to the names of each SubHeader
	void changeKey(std::string oldKey, std::string newKey);	///< Change one of the keys in the map
	void addPair(std::string key, std::string val);			///< Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int t);	///< Adds a pair object and asserts a type
	void setName(std::string name);							///< Set the name of the Header
	void setAlias(std::string alias);						///< Set the alias of the header, if any
	void setNameAliasPair(std::string n, std::string a, int s);	///< Set the name, alias, and state for the header
	void setState(int state);								///< Set the state of the header, if any
	void DisplayContents();									///< Display the contents of the header object
	
	void addSubKey(std::string key);						///< Adds a key to the SubHeader Map
	void copyAnchor2Alias(std::string alias, SubHeader &ref);///< Find the anchor in the map, and copy to the Header reference given
	
	int size();												///< Return the size of the Sub_Map
	std::string getName();									///< Return the name of the header
	std::string getAlias();									///< Return the alias of the header
	int getState();											///< Return the state of the header
	bool isAlias();											///< Returns true if the header is an alias
	bool isAnchor();										///< Returns true if the header is an anchor
	SubHeader& getAnchoredSub(std::string alias);			///< Returns reference to the anchored subheader, if any
	
protected:
	
private:
	std::map<std::string, SubHeader> Sub_Map;			///< Map of the contained subheaders in the main header
	
};

/// Object for the various documents in the yaml file
/** C++ Object for the documents in a yaml input file as denoted by a Key: followed by --- (three dashes) and ending
	with a ... (three dots). A single yaml file can have multiple document structures and each document structure can
	have multiple headers (which have sub-headers and key-values) and key-value-pairs. This is the larges single
	object in the yaml file itself. 
 
	Just like Header, this object also inherits from SubHeader and therefore has access to its protected members. 
	You can use access to those members to establish the KeyValuePairs in the Document, name the Document, and
	give the Document an alias or anchor value. */
class Document : SubHeader
{
public:
	Document();											///< Default constructor
	~Document();										///< Default destructor
	Document(const Document &doc);						///< Copy constructor
	Document(std::string name);							///< Constructor by name
	Document(const KeyValueMap &map);					///< Constructor by existing map
	Document(std::string name, const KeyValueMap &map);	///< Constructor by name and map
	Document(std::string key, const Header &head);		///< Constructor by single header
	
	Document& operator=(const Document &doc);				///< Equals overload
	ValueTypePair& operator[](const std::string key);		///< Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	///< Return the ValueType at the given key
	Header& operator()(const std::string key);				///< Return the Header reference at the given key
	Header operator()(const std::string key) const;			///< Return the Header at the given key
	
	std::map<std::string, Header> & getHeadMap();			///< Return the reference to the Header Map
	KeyValueMap & getDataMap();								///< Return the reference to the KeyValueMap
	Header & getHeader(std::string key);					///< Return reference to the Header in map at the key
	std::map<std::string, Header>::const_iterator end() const;	///< Returns a const iterator pointing to the end of the list
	std::map<std::string, Header>::iterator end();				///< Returns an iterator pointing to the end of the list
	std::map<std::string, Header>::const_iterator begin() const;///< Returns a const iterator pointing to the begining of the list
	std::map<std::string, Header>::iterator begin();			///< Returns an iterator pointing to the begining of the list
	
	void clear();											///< Clear out info in the Document
	void resetKeys();										///< Set all keys in the map to match names of the headers
	void changeKey(std::string oldKey, std::string newKey);	///< Change a given oldKey in the header map to the newKey given
	void revalidateAllKeys();								///< Resets and validates keys in header and subheader maps
	void addPair(std::string key, std::string val);			///< Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int t);	///< Adds a pair object and asserts a type
	void setName(std::string name);							///< Set the name of the Document
	void setAlias(std::string alias);						///< Set the alias of the Document
	void setNameAliasPair(std::string n, std::string a, int s);///< Set the name, alias, and state of the document
	void setState(int state);								///< Set the state of the Document
	void DisplayContents();									///< Display the contents of the Document
	
	void addHeadKey(std::string key);						///< Add a key to the Header without a header object
	void copyAnchor2Alias(std::string alias, Header &ref);	///< Find the anchor in the map, and copy to the Header reference given
	
	int size();												///< Return the size of the header map
	std::string getName();									///< Return the name of the document
	std::string getAlias();									///< Return the alias of the document
	int getState();											///< Return the state of the document
	bool isAlias();											///< Returns true if the document is an alias
	bool isAnchor();										///< Returns true if the document is an anchor
	Header& getAnchoredHeader(std::string alias);			///< Returns reference to the anchored header, if any
	Header& getHeadFromSubAlias(std::string alias);			///< Returns reference to the Header that contains a Sub with the given alias
	
protected:
	
private:
	std::map<std::string, Header> Head_Map;				///< Map of headers contained within the document
	
};

/// Object for the entire yaml file holding all documents, header, sub-headers, keys, and values
/** C++ Object for the yaml file. This object holds a map of all Documents in the yaml file. Each document
	holds a map of Key-values and Headers. The Headers hold maps of Key-values and SubHeaders, and each
	SubHeader can hold more Key-values. 
 
	This object is used to represent a digital and queryable structure for all information contained within
	a yaml file. There are some limitations to what can be held here, and those limitations are based on the
	limitations in the C-YAML Library token parser. The main limitation is that the deepest level of allowable
	recursion in the file is SubHeader. Meaning that you are not allowed to have Sub-SubHeaders underneath a
	SubHeader object. This imposes a hard limit to number of nested lists that can be in a single Document
	object. 
 
	When using yaml_cpp_class, this object will generally be what you query into to get the information from
	your yaml input files. From this object, functions and operators are provided to give you the capability
	of querying down into any allowable level of the file by the keys that are were used in the file. Be sure 
	that you are querying the correct objects by the correct keys, otherwise errors and exceptions will be
	thrown. */
class YamlWrapper
{
public:
	YamlWrapper();										///< Default constructor
	~YamlWrapper();										///< Default destructor
	YamlWrapper(const YamlWrapper &yaml);				///< Copy constructor
	YamlWrapper(std::string key, const Document &doc);	///< Constructor by a single document
	
	YamlWrapper& operator=(const YamlWrapper &yaml);	///< Equals overload
	Document& operator()(const std::string key);		///< Return the Document reference at the given key
	Document operator()(const std::string key) const;	///< Return the Document at the given key
	
	std::map<std::string, Document> & getDocMap();		///< Return reference to the document map
	Document& getDocument(std::string key);				///< Return reference to the document at the key
	std::map<std::string, Document>::const_iterator end() const;	///< Returns a const iterator pointing to the end of the list
	std::map<std::string, Document>::iterator end();				///< Returns an iterator pointing to the end of the list
	std::map<std::string, Document>::const_iterator begin() const;	///< Returns a const iterator pointing to the begining of the list
	std::map<std::string, Document>::iterator begin();				///< Returns an iterator pointing to the begining of the list
	
	void clear();										///< Clear out the yaml object
	void resetKeys();									///< Resets all the keys in DocumentMap to match document names
	void changeKey(std::string oldKey, std::string newKey);	///< Change a given oldKey in the map to the newKey given
	void revalidateAllKeys();							///< Resets and validates all keys in the structure
	void DisplayContents();								///< Display the contents of the wrapper
	
	void addDocKey(std::string key);					///< Add a key to the document map
	void copyAnchor2Alias(std::string alias, Document &ref);///< Find the anchor in the map, and copy to the Document reference given
	
	int size();											///< Return the size of the document map
	Document& getAnchoredDoc(std::string alias);		///< Return the reference to the document that is anchored with the given alias
	Document& getDocFromHeadAlias(std::string alias);	///< Return reference to the document that contains the header with the given alias
	Document& getDocFromSubAlias(std::string alias);	///< Return reference to the document that contains the subheader with the given alias
	
protected:
	
private:
	std::map<std::string, Document> Doc_Map;			///< Map of the documents contained within the wrapper
	
};

/// Primary object used when reading and digitally storing yaml files
/** C++ Object that holds the YamlWrapper object and the C-YAML objects necessary for reading and parsing a yaml
	formatted file. This is the primary object that users are expected to work with when using yaml_wrapper.h to
	read input files. It contains functions necessary to setup a read instance, read and parse the input, place
	the parsed input results into the digital YamlWrapper object, and allow the user to query that object. 
 
	The two main functions that the typical user will need are: (i) executeYamlRead() and (ii) getYamlWrapper. 
	Make sure that the read function was called prior to querying the YamlWrapper structure. Do not call the
	cleanup() function if using executeYamlRead(). That function will be automattically called after the 
	read is complete. */
class yaml_cpp_class
{
public:
	yaml_cpp_class();									///< Default constructor
	~yaml_cpp_class();									///< Default destructor
	
	int setInputFile(const char *file);					///< Set the input file to be read
	int readInputFile();								///< Reads through input file and stores into YamlWrapper
	int cleanup();										///< Deletes yaml_c objects and closes the input file
	int executeYamlRead(const char *file);				///< Runs the full execution of initialization, reading, and cleaning
	YamlWrapper& getYamlWrapper();						///< Returns reference to the YamlWrapper Object
	
	void DisplayContents();								///< Print out the contents of the read to the console window
	
protected:
	
private:
	YamlWrapper yaml_wrapper;							///< YamlWrapper object where digital file is stored
	FILE *input_file;									///< Function pointer to the yaml formatted file
	const char *file_name;								///< Name of the file to be parsed and read
	yaml_parser_t token_parser;							///< C-YAML parser object for token based parsing
	yaml_token_t current_token;							///< C-YAML token object for the current token in the file
	yaml_token_t previous_token;						///< C-YAML token object for the previous token in the file
	
};

/// Function to return an all lower case string based on the passed argument
/** This function will copy the input paramter and convert that copy to
	all lower case. The copy is then returned and can be checked against valid
	or allowed strings.
 
 \param input string to copy and convert to lower case */
std::string allLower(const std::string &input);

/// Function to return true if the given argument is an even number
bool isEven(int n);

/// Function to run tests on all the objects that yaml_cpp_class
/** This test function is currently NOT callable from the UI.*/
int YAML_WRAPPER_TESTS();

/// Function to run a test read for the yaml_cpp_class on a given file
/** This test/executable function is currently NOT callable from the UI.*/
int YAML_CPP_TEST(const char *file);

#endif
