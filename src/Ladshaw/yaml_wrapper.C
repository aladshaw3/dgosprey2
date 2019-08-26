/*!
 *  \file yaml_wrapper.cpp yaml_wrapper.h
 *	\brief C++ Wrapper for the C-YAML Library
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


#include "yaml_wrapper.h"

/*
  =========================================== ValueTypePair Class ==============================================
 */

//Default constructor
ValueTypePair::ValueTypePair()
{
	type = UNKNOWN;
	Value_Type = std::make_pair(" ", type);
}

//Default Destructor
ValueTypePair::~ValueTypePair()
{
	Value_Type.first.clear();
}

//Construction by pair
ValueTypePair::ValueTypePair(const std::pair<std::string, int> &vt)
{
	Value_Type = vt;
	int t = vt.second;
	if (t != STRING && t != BOOLEAN && t != DOUBLE && t != INT && t != UNKNOWN)
	{
		mError(invalid_type);
		t = UNKNOWN;
	}
	type = t;
}

//Construction by string value and int type
ValueTypePair::ValueTypePair(std::string value, int t)
{
	Value_Type = std::make_pair(value, t);
	if (t != STRING && t != BOOLEAN && t != DOUBLE && t != INT && t != UNKNOWN)
	{
		mError(invalid_type);
		t = UNKNOWN;
	}
	type = t;
}

//Copy constructor for ValueTypePair
ValueTypePair::ValueTypePair(const ValueTypePair &vt)
{
	Value_Type = vt.Value_Type;
	type = vt.type;
}

//Equals overload
ValueTypePair& ValueTypePair::operator=(const ValueTypePair &vt)
{
	if (this == &vt)
		return *this;
	else
	{
		this->Value_Type = vt.Value_Type;
		this->type = vt.type;
		return *this;
	}
}

//Edit the value to the pair
void ValueTypePair::editValue(std::string value)
{
	this->Value_Type.first = value;
}

//Edit the pair
void ValueTypePair::editPair(std::string value, int t)
{
	this->Value_Type.first = value;
	if (t != STRING && t != BOOLEAN && t != DOUBLE && t != INT && t != UNKNOWN)
	{
		mError(invalid_type);
		t = UNKNOWN;
	}
	this->Value_Type.second = t;
	this->type = t;
}

//Find out the type of data
void ValueTypePair::findType()
{
	//Check to see if we have a number
	bool haveDOT = false;
	bool haveE = false;
	bool havePlusMinus = false;
	bool isNumber = true;
	
	//Check to see if it is a bool
	std::string copy = this->Value_Type.first;
	
	//Loop through the string to determine whether or not this is a number
	for (int i=0; i<this->Value_Type.first.length(); i++)
	{
		if (!isdigit(this->Value_Type.first[i]))
		{
			if ( !(this->Value_Type.first[i] == '-' || this->Value_Type.first[i] == '+') && i==0 )
			{
				isNumber = false;
				break;
			}
			else if ( this->Value_Type.first[i] == '.')
			{
				if (haveDOT == true)
				{
					isNumber = false;
					break;
				}
				else
					haveDOT = true;
			}
			else if ( (this->Value_Type.first[i] == 'E' || this->Value_Type.first[i] == 'e') && i>0  && this->Value_Type.first[i-1] != ' ')
			{
				if (haveE == true)
				{
					isNumber = false;
					break;
				}
				else
				{
					haveE = true;
					if ( i == this->Value_Type.first.length()-1)
					{
						isNumber = false;
						break;
					}
					else
					{
						i++;
						if (!isdigit(this->Value_Type.first[i]))
						{
							if ( this->Value_Type.first[i] == '-' || this->Value_Type.first[i] == '+')
							{
								if (havePlusMinus == true)
								{
									isNumber = false;
									break;
								}
								else
									havePlusMinus = true;
							}
							else
							{
								isNumber = false;
								break;
							}
						}
						else
						{
							//continue
						}
					}
				}
			}
			else if (i!=0 && (this->Value_Type.first[i] == '-' || this->Value_Type.first[i] == '+'))
			{
				isNumber = false;
				break;
			}
			else
			{
				//continue
			}
		}
		else
		{
			//continue
		}
	} //END LOOP
	
	//If it is a number, set the type and return
	if (isNumber == true && haveDOT == false && haveE == false)
	{
		this->type = INT;
	}
	else if (isNumber == true && (haveDOT == true || haveE == true))
	{
		this->type = DOUBLE;
	}
	else
	{
		//Determine whether we have a boolean or a word
		for (int i=0; i<this->Value_Type.first.length(); i++)
		{
			copy[i] = tolower(copy[i]);
		}
		if (copy == "false" || copy == "true")
			this->type = BOOLEAN;
		else
			this->type = STRING;
	}
	this->Value_Type.second = this->type;
}

//Force a specific type
void ValueTypePair::assertType(int t)
{
	
	if (!(t==UNKNOWN || t==STRING || t==BOOLEAN || t==DOUBLE || t==INT))
	{
		mError(invalid_type);
		return;
	}
	this->type = t;
	this->Value_Type.second = t;
}

//Display the pair
void ValueTypePair::DisplayPair()
{
	std::cout << "Value = " << this->Value_Type.first << "\tType = ";
	if (this->type == UNKNOWN)
		std::cout << "UNKNOWN\n";
	else if (this->type == STRING)
		std::cout << "STRING\n";
	else if (this->type == BOOLEAN)
		std::cout << "BOOLEAN\n";
	else if (this->type == DOUBLE)
		std::cout << "DOUBLE\n";
	else if (this->type == INT)
		std::cout << "INT\n";
	else
		mError(invalid_type);
	std::cout << std::endl;
}

//Return the string of value in pair
std::string ValueTypePair::getString()
{
	if (this->Value_Type.second != STRING && this->Value_Type.second != UNKNOWN)
		mError(invalid_type);
	return this->Value_Type.first;
}

//Returns the boolean of the value in pair
bool ValueTypePair::getBool()
{
	if (this->Value_Type.second != BOOLEAN)
	{
		mError(invalid_type);
		return false;
	}
	else
	{
		std::string copy = this->Value_Type.first;
		//Determine whether we have a boolean or a word
		for (int i=0; i<this->Value_Type.first.length(); i++)
		{
			copy[i] = tolower(copy[i]);
		}
		if (copy == "false")
			return false;
		else
			return true;
	}
}

//Returns the double of the value
double ValueTypePair::getDouble()
{
	if (this->Value_Type.second != DOUBLE && this->Value_Type.second != INT)
		mError(invalid_type);
	return atof(this->Value_Type.first.c_str());
}

//Returns the int of the value
int ValueTypePair::getInt()
{
	if (this->Value_Type.second != INT && this->Value_Type.second != DOUBLE)
		mError(invalid_type);
	return atoi(this->Value_Type.first.c_str());
}

//Return the value of the pair as a string
std::string ValueTypePair::getValue()
{
	return this->Value_Type.first;
}

//Return the type of the pair
int ValueTypePair::getType()
{
	return this->type;
}

//Return the actual pair
std::pair<std::string,int>& ValueTypePair::getPair()
{
	return this->Value_Type;
}

/*
 =========================================== END ValueTypePair Class ==============================================
 */

/*
 =========================================== KeyValueMap Class ==============================================
 */

//Default constructor
KeyValueMap::KeyValueMap()
{
}

//Default destructor
KeyValueMap::~KeyValueMap()
{
	Key_Value.clear();
}

//Construct from a string map
KeyValueMap::KeyValueMap(const std::map<std::string, std::string> &map)
{
	for (auto &x: map)
	{
		ValueTypePair temp(x.second,UNKNOWN);
		Key_Value.insert( Key_Value.end(), std::pair<std::string,ValueTypePair> (x.first,temp) );
	}
}

//Construct with one element
KeyValueMap::KeyValueMap(std::string key, std::string value)
{
	ValueTypePair temp(value,UNKNOWN);
	Key_Value.insert(Key_Value.end(), std::pair<std::string,ValueTypePair>(key,temp) );
}

//Copy constructor
KeyValueMap::KeyValueMap(const KeyValueMap &map)
{
	Key_Value = map.Key_Value;
}

//Equals overload for KeyValueMap
KeyValueMap& KeyValueMap::operator=(const KeyValueMap &map)
{
	if (this == &map)
		return *this;
	else
	{
		this->Key_Value.clear();
		this->Key_Value = map.Key_Value;
		return *this;
	}
}

//Reference fetch overload
ValueTypePair& KeyValueMap::operator[](const std::string key)
{
	try
	{
		return this->Key_Value.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::iterator it;
		it = this->Key_Value.end();
		return it->second;
	}
}

//Access operator overload
ValueTypePair KeyValueMap::operator[](const std::string key) const
{
	try
	{
		return this->Key_Value.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::const_iterator it;
		it = this->Key_Value.end();
		return it->second;
	}
}

//Returns a reference to the map
std::map<std::string, ValueTypePair >& KeyValueMap::getMap()
{
	return this->Key_Value;
}

//Returns a const iterator pointing at the end of the list
std::map<std::string, ValueTypePair>::const_iterator KeyValueMap::end() const
{
	std::map<std::string, ValueTypePair>::const_iterator it;
	it = this->Key_Value.end();
	return it;
}

//Returns a iterator pointing at the end of the list
std::map<std::string, ValueTypePair>::iterator KeyValueMap::end()
{
	std::map<std::string, ValueTypePair>::iterator it;
	it = this->Key_Value.end();
	return it;
}

//Returns a const iterator pointing at the beginning of the list
std::map<std::string, ValueTypePair>::const_iterator KeyValueMap::begin() const
{
	std::map<std::string, ValueTypePair>::const_iterator it;
	it = this->Key_Value.begin();
	return it;
}

//Returns a iterator pointing at the beginning of the list
std::map<std::string, ValueTypePair>::iterator KeyValueMap::begin()
{
	std::map<std::string, ValueTypePair>::iterator it;
	it = this->Key_Value.begin();
	return it;
}

//Clear the map object
void KeyValueMap::clear()
{
	this->Key_Value.clear();
}

//Add a black key to the map
void KeyValueMap::addKey(std::string key)
{
	this->addPair(key, " ", UNKNOWN);
}

//Edit a value for a pre-existing key
void KeyValueMap::editValue4Key(std::string val, std::string key)
{
	try
	{
		this->Key_Value.at(key).editValue(val);
		this->Key_Value.at(key).findType();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
	}
}

//Edit a value for a key and assert the type
void KeyValueMap::editValue4Key(std::string val, int t,std::string key)
{
	try
	{
		this->Key_Value.at(key).editPair(val, t);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
	}
}

//Add a paired object to the map
void KeyValueMap::addPair(std::string key, ValueTypePair val)
{
	this->Key_Value.insert(this->Key_Value.end(), std::pair<std::string,ValueTypePair>(key, val) );
}

//Add a paired object with only strings
void KeyValueMap::addPair(std::string key, std::string val)
{
	ValueTypePair temp(val,UNKNOWN);
	temp.findType();
	this->Key_Value.insert(this->Key_Value.end(), std::pair<std::string,ValueTypePair>(key, temp) );
}

//Add a paired object and assert its type
void KeyValueMap::addPair(std::string key, std::string val, int t)
{
	ValueTypePair temp(val,t);
	this->Key_Value.insert(this->Key_Value.end(), std::pair<std::string,ValueTypePair>(key, temp) );
}

//Find the data type at the key
void KeyValueMap::findType(std::string key)
{
	try
	{
		this->Key_Value.at(key).findType();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
	}
}

//Assert a data type by the key
void KeyValueMap::assertType(std::string key, int t)
{
	try
	{
		this->Key_Value.at(key).assertType(t);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
	}
}

//Find all types for all data in the map
void KeyValueMap::findAllTypes()
{
	for (auto &x: this->Key_Value)
	{
		x.second.findType();
	}
}

//Display the contents of the map
void KeyValueMap::DisplayMap()
{
	std::cout << "KeyValueMap Contents\n---------------------\n";
	for (auto &x: this->Key_Value)
	{
		std::cout << "Key: " << x.first << "\tValue: " << x.second.getValue() << "\tType: ";
		if (x.second.getType() == UNKNOWN)
			std::cout << "UNKNOWN" << std::endl;
		else if (x.second.getType() == STRING)
			std::cout << "STRING" << std::endl;
		else if (x.second.getType() == BOOLEAN)
			std::cout << "BOOLEAN" << std::endl;
		else if (x.second.getType() == DOUBLE)
			std::cout << "DOUBLE" << std::endl;
		else if (x.second.getType() == INT)
			std::cout << "INT" << std::endl;
		else
			mError(invalid_type);
	}
	std::cout << std::endl;
}

//Return the size of the map object
int KeyValueMap::size()
{
	return (int) this->Key_Value.size();
}

//Get the string for the key
std::string KeyValueMap::getString(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getString();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return "NO_KEY";
	}
}

//Get the bool for the key
bool KeyValueMap::getBool(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getBool();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return false;
	}
}

//Get the double for the key
double KeyValueMap::getDouble(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getDouble();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return 0;
	}
}

//Get the int for the key
int KeyValueMap::getInt(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getInt();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return 0;
	}
}

//Get the value of the key
std::string KeyValueMap::getValue(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getValue();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return "NO_KEY";
	}
}

//Get the type for the key
int KeyValueMap::getType(std::string key)
{
	try
	{
		return this->Key_Value.at(key).getType();
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		return UNKNOWN;
	}
}

//Get the pair at the key
ValueTypePair& KeyValueMap::getPair(std::string key)
{
	try
	{
		return this->Key_Value.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::iterator it;
		it = this->Key_Value.end();
		return it->second;
	}
}

/*
 =========================================== END KeyValueMap Class ==============================================
 */

/*
 =========================================== SubHeader Class ==============================================
 */

//Default constructor for SubHeader
SubHeader::SubHeader()
{
	state = NONE;
	alias = " ";
	name = " ";
}

//Default destructor for SubHeader
SubHeader::~SubHeader()
{
	Data_Map.clear();
	alias.clear();
	name.clear();
}

//Copy constructor
SubHeader::SubHeader(const SubHeader &subheader)
{
	Data_Map = subheader.Data_Map;
	alias = subheader.alias;
	name = subheader.name;
	state = subheader.state;
}

//Construction be existing map
SubHeader::SubHeader(const KeyValueMap &map)
{
	Data_Map = map;
	alias = " ";
	name = " ";
	state = NONE;
}

//Construction by name
SubHeader::SubHeader(std::string n)
{
	alias = " ";
	name = n;
	state = NONE;
}

//Construction by name and map
SubHeader::SubHeader(std::string n, const KeyValueMap &map)
{
	Data_Map = map;
	alias = " ";
	name = n;
	state = NONE;
}

//Equals overload
SubHeader& SubHeader::operator=(const SubHeader &sub)
{
	if (this == &sub)
		return *this;
	else
	{
		this->Data_Map.clear();
		this->Data_Map = sub.Data_Map;
		this->alias = sub.alias;
		this->name = sub.name;
		this->state = sub.state;
		return *this;
	}
}

//Reference fetch overload
ValueTypePair& SubHeader::operator[](const std::string key)
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Access operator overload
ValueTypePair SubHeader::operator[](const std::string key) const
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::const_iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Return reference to the KeyValue Map
KeyValueMap& SubHeader::getMap()
{
	return this->Data_Map;
}

//Clear out the data structure
void SubHeader::clear()
{
	this->Data_Map.clear();
	this->name.clear();
	this->alias.clear();
	this->state = NONE;
}

//Adds a pair to the KeyValueMap in the subheader
void SubHeader::addPair(std::string key, std::string val)
{
	this->Data_Map.addPair(key, val);
}

//Adds a pair to the KeyValueMap and asserts type
void SubHeader::addPair(std::string key, std::string val, int type)
{
	this->Data_Map.addPair(key, val, type);
}

//Set the name of the subheader
void SubHeader::setName(std::string n)
{
	this->name = n;
}

//Set the alias without a type
void SubHeader::setAlias(std::string alias)
{
	this->alias = alias;
	this->state = NONE;
}

//Set the alais of the subheader
void SubHeader::setAlias(std::string a, int s)
{
	this->alias = a;
	if (s != ANCHOR && s != ALIAS && s != NONE)
	{
		mError(invalid_type);
		s = NONE;
	}
	this->state = s;
}

//Set both name and alias
void SubHeader::setNameAliasPair(std::string n, std::string a, int s)
{
	this->setName(n);
	this->setAlias(a,s);
}

//Set the state
void SubHeader::setState(int s)
{
	if (s != ANCHOR && s != ALIAS && s != NONE)
	{
		mError(invalid_type);
		s = NONE;
	}
	this->state = s;
}

//Display the contents of the subheader
void SubHeader::DisplayContents()
{
	std::cout << "SubHeader Class\nName: " << this->getName();
	if (this->isAnchor())
		std::cout << " -> ANCHOR: " << this->getAlias() << "\n----------------------\n";
	else if (this->isAlias())
		std::cout << " -> ALIAS: " << this->getAlias() << "\n----------------------\n";
	else
		std::cout << "\n----------------------\n";
	this->Data_Map.DisplayMap();
}

//Return the name of the subheader
std::string SubHeader::getName()
{
	return this->name;
}

//Return the alias of the subheader
std::string SubHeader::getAlias()
{
	return this->alias;
}

//Return true if subheader is an alias
bool SubHeader::isAlias()
{
	if (this->state == ALIAS)
		return true;
	else
		return false;
}

//Return true if subheader is an anchor
bool SubHeader::isAnchor()
{
	if (this->state == ANCHOR)
		return true;
	else
		return false;
}

//Return the state of the subheader
int SubHeader::getState()
{
	return this->state;
}

/*
 =========================================== END SubHeader Class ==============================================
 */

/*
 =========================================== Header Class ==============================================
 */

//Default constructor
Header::Header()
{
	state = NONE;
	alias = " ";
	name = " ";
}

//Default destructor
Header::~Header()
{
	Data_Map.clear();
	name.clear();
	alias.clear();
	Sub_Map.clear();
}

//Copy Constructor
Header::Header(const Header &head)
{
	Data_Map = head.Data_Map;
	Sub_Map = head.Sub_Map;
	name = head.name;
	alias = head.alias;
	state = head.state;
}

//Constructor by name
Header::Header(std::string n)
{
	name = n;
	alias = " ";
	state = NONE;
}

//Constructor by map
Header::Header(const KeyValueMap &map)
{
	Data_Map = map;
	name = " ";
	alias = " ";
	state = NONE;
}

//Constructor by name and map
Header::Header(std::string n, const KeyValueMap &map)
{
	Data_Map = map;
	name = n;
	alias = " ";
	state = NONE;
}

//Constructor by single subheader object
Header::Header(std::string key, const SubHeader &sub)
{
	name = " ";
	alias = " ";
	state = NONE;
	Sub_Map.insert(Sub_Map.end(), std::pair<std::string, SubHeader>(key, sub) );
}

//Equals overload
Header& Header::operator=(const Header &head)
{
	if (this == &head)
	{
		return *this;
	}
	else
	{
		this->Data_Map.clear();
		this->Sub_Map.clear();
		this->Data_Map = head.Data_Map;
		this->Sub_Map = head.Sub_Map;
		this->name = head.name;
		this->alias = head.alias;
		this->state	= head.state;
		return *this;
	}
}

//Reference fetch overload
ValueTypePair& Header::operator[](const std::string key)
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Access operator overload
ValueTypePair Header::operator[](const std::string key) const
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::const_iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Reference fetch overload
SubHeader& Header::operator()(const std::string key)
{
	try
	{
		return this->Sub_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, SubHeader>::iterator it;
		it = this->Sub_Map.end();
		return it->second;
	}
}

//Access operator overload
SubHeader Header::operator()(const std::string key) const
{
	try
	{
		return this->Sub_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, SubHeader>::const_iterator it;
		it = this->Sub_Map.end();
		return it->second;
	}
}

//Return reference to the SubMap
std::map<std::string, SubHeader> & Header::getSubMap()
{
	return this->Sub_Map;
}

//Return reference to the DataMap
KeyValueMap & Header::getDataMap()
{
	return this->getMap();
}

//Return reference to the SubHeader object at the given key
SubHeader & Header::getSubHeader(std::string key)
{
	try
	{
		return this->Sub_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, SubHeader>::iterator it;
		it = this->Sub_Map.end();
		return it->second;
	}
}

//Return a constant iterator for the end of the custom map
std::map<std::string, SubHeader>::const_iterator Header::end() const
{
	std::map<std::string, SubHeader>::const_iterator it;
	it = this->Sub_Map.end();
	return it;
}

//Return an iterator for the end of the custom map
std::map<std::string, SubHeader>::iterator Header::end()
{
	std::map<std::string, SubHeader>::iterator it;
	it = this->Sub_Map.end();
	return it;
}

//Return a constant iterator for the beginning of the custom map
std::map<std::string, SubHeader>::const_iterator Header::begin() const
{
	std::map<std::string, SubHeader>::const_iterator it;
	it = this->Sub_Map.begin();
	return it;
}

//Return an iterator for the beginning of the custom map
std::map<std::string, SubHeader>::iterator Header::begin()
{
	std::map<std::string, SubHeader>::iterator it;
	it = this->Sub_Map.begin();
	return it;
}

//Clear out the memory
void Header::clear()
{
	this->Sub_Map.clear();
	this->Data_Map.clear();
	this->name.clear();
	this->alias.clear();
	this->state = NONE;
}

//Reset the keys of the SubMap to reflect the names of the SubHeaders
void Header::resetKeys()
{
	std::string newKey;
	std::string oldKey;
	for (auto &x: this->Sub_Map)
	{
		oldKey = x.first;
		newKey = x.second.getName();
		std::swap(this->Sub_Map[newKey],x.second);
		this->Sub_Map.erase(x.first);
	}
}

//Change one of the keys in the map
void Header::changeKey(std::string oldKey, std::string newKey)
{
	std::map<std::string, SubHeader>::const_iterator it = this->Sub_Map.find(oldKey);
	if (it != this->Sub_Map.end())
	{
		SubHeader temp(it->second);
		this->Sub_Map.erase(it);
		this->Sub_Map.insert(this->Sub_Map.end(), std::pair<std::string, SubHeader>(newKey,temp) );
	}
	else
	{
		mError(key_not_found);
	}
}

//Add a pair to the KeyValueMap of the header
void Header::addPair(std::string key, std::string val)
{
	this->Data_Map.addPair(key, val);
}

//Add a pair to the KeyValueMap of the header and assert type
void Header::addPair(std::string key, std::string val, int t)
{
	this->Data_Map.addPair(key, val, t);
}

//Set the name of the header
void Header::setName(std::string n)
{
	this->name = n;
}

//Set the alias of the header
void Header::setAlias(std::string a)
{
	this->alias = a;
}

//Set the name alias and state of the header
void Header::setNameAliasPair(std::string n, std::string a, int s)
{
	this->SubHeader::setNameAliasPair(n, a, s);
}

//Set the state of the header
void Header::setState(int s)
{
	this->SubHeader::setState(s);
}

//Display the contents of the Header object
void Header::DisplayContents()
{
	std::cout << "Header Class\nName: " << this->getName();
	if (this->isAnchor())
		std::cout << " -> ANCHOR: " << this->getAlias() << "\n----------------------\n";
	else if (this->isAlias())
		std::cout << " -> ALIAS: " << this->getAlias() << "\n----------------------\n";
	else
		std::cout << "\n----------------------\n";
	std::cout << "Contents of this header include " << this->size() << " SubHeaders and " << this->Data_Map.size() << " KeyValue Pairs...\n\n";
	
	//Display contents of KeyValueMap
	if (this->Data_Map.size() > 0)
		this->Data_Map.DisplayMap();
	
	//Loop through the SubHeaders
	for (auto &x: this->Sub_Map)
	{
		x.second.DisplayContents();
	}
}

//Add a SubHeader key to the Header
void Header::addSubKey(std::string key)
{
	SubHeader temp(key);
	this->Sub_Map.insert(this->Sub_Map.end(), std::pair<std::string, SubHeader>(key,temp) );
}

//Find and copy the contents of the SubHeader Anchor to the Alias SubHeader
void Header::copyAnchor2Alias(std::string alias, SubHeader &ref)
{
	std::string key = ref.getName();
	ref = getAnchoredSub(alias);
	ref.setNameAliasPair(key, alias, ALIAS);
	
}


//Return the size of the Sub_Map
int Header::size()
{
	return (int) this->Sub_Map.size();
}

//Return the name of the header
std::string Header::getName()
{
	return this->name;
}

//Return the alias of the header
std::string Header::getAlias()
{
	return this->alias;
}

//Return the state of the header
int Header::getState()
{
	return this->state;
}

//Return true if header is an alias
bool Header::isAlias()
{
	if (this->getState() == ALIAS)
		return true;
	else
		return false;
}

//Return true if header is an anchor
bool Header::isAnchor()
{
	if (this->getState() == ANCHOR)
		return true;
	else
		return false;
}

//Return reference to the anchored subheader
SubHeader& Header::getAnchoredSub(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Sub_Map)
	{
		if (x.second.isAnchor())
		{
			if (alias == x.second.getAlias())
			{
				foundit = true;
				current = x.first;
				break;
			}
		}
		else
		{
			current = x.first;
		}
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getSubHeader(current);
}


/*
 =========================================== END Header Class ==============================================
 */

/*
 =========================================== Document Class ==============================================
 */

//Default constructor
Document::Document()
{
	name = " ";
	alias = " ";
	state = NONE;
}

//Default destructor
Document::~Document()
{
	name.clear();
	alias.clear();
	Head_Map.clear();
	Data_Map.clear();
}

//Copy constructor
Document::Document(const Document &doc)
{
	name = doc.name;
	alias = doc.alias;
	state = doc.state;
	Head_Map = doc.Head_Map;
	Data_Map = doc.Data_Map;
}

//Constructor by name
Document::Document(std::string n)
{
	name = n;
	alias = " ";
	state = NONE;
}

//Constructor by existing map
Document::Document(const KeyValueMap & map)
{
	name = " ";
	alias = " ";
	state = NONE;
	Data_Map = map;
}

//Constructor by name and map
Document::Document(std::string n, const KeyValueMap & map)
{
	name = n;
	alias = " ";
	state = NONE;
	Data_Map = map;
}

//Constructor by a key and header
Document::Document(std::string key, const Header & head)
{
	name = " ";
	alias = " ";
	state = NONE;
	Head_Map.insert(Head_Map.end(), std::pair<std::string, Header>(key, head) );
}

//Equals overload
Document& Document::operator=(const Document &doc)
{
	if (this == &doc)
		return *this;
	else
	{
		this->Data_Map.clear();
		this->Head_Map.clear();
		this->Data_Map = doc.Data_Map;
		this->Head_Map = doc.Head_Map;
		this->name = doc.name;
		this->alias = doc.alias;
		this->state = doc.state;
		return *this;
	}
}

//Reference fetch overload
ValueTypePair& Document::operator[](const std::string key)
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Access operator overload
ValueTypePair Document::operator[](const std::string key) const
{
	try
	{
		return this->Data_Map[key];
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, ValueTypePair>::const_iterator it;
		it = this->Data_Map.end();
		return it->second;
	}
}

//Reference fetch overload
Header& Document::operator()(const std::string key)
{
	try
	{
		return this->Head_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Header>::iterator it;
		it = this->Head_Map.end();
		return it->second;
	}
}

//Access operator overload
Header Document::operator()(const std::string key) const
{
	try
	{
		return this->Head_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Header>::const_iterator it;
		it = this->Head_Map.end();
		return it->second;
	}
}

//Return reference to the header map
std::map<std::string, Header> & Document::getHeadMap()
{
	return this->Head_Map;
}

//Return reference to the key-value map
KeyValueMap & Document::getDataMap()
{
	return this->Data_Map;
}

//Return reference to the header in the map at the given key
Header & Document::getHeader(std::string key)
{
	try
	{
		return this->Head_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Header>::iterator it;
		it = this->Head_Map.end();
		return it->second;
	}
}

//Return a constant iterator for the end of the custom map
std::map<std::string, Header>::const_iterator Document::end() const
{
	std::map<std::string, Header>::const_iterator it;
	it = this->Head_Map.end();
	return it;
}

//Return an iterator for the end of the custom map
std::map<std::string, Header>::iterator Document::end()
{
	std::map<std::string, Header>::iterator it;
	it = this->Head_Map.end();
	return it;
}

//Return a constant iterator for the beginning of the custom map
std::map<std::string, Header>::const_iterator Document::begin() const
{
	std::map<std::string, Header>::const_iterator it;
	it = this->Head_Map.begin();
	return it;
}

//Return an iterator for the beginning of the custom map
std::map<std::string, Header>::iterator Document::begin()
{
	std::map<std::string, Header>::iterator it;
	it = this->Head_Map.begin();
	return it;
}

//Clear out the memory in the Document
void Document::clear()
{
	this->name.clear();
	this->alias.clear();
	this->state = NONE;
	this->Data_Map.clear();
	this->Head_Map.clear();
}

//Reset all keys in the map to match names given to the headers
void Document::resetKeys()
{
	std::string newKey;
	std::string oldKey;
	for (auto &x: this->Head_Map)
	{
		oldKey = x.first;
		newKey = x.second.getName();
		std::swap(this->Head_Map[newKey],x.second);
		this->Head_Map.erase(x.first);
	}
}

//Change the given oldKey to the newKey in the map, if oldKey exists
void Document::changeKey(std::string oldKey, std::string newKey)
{
	std::map<std::string, Header>::const_iterator it = this->Head_Map.find(oldKey);
	if (it != this->Head_Map.end())
	{
		Header temp(it->second);
		this->Head_Map.erase(it);
		this->Head_Map.insert(this->Head_Map.end(), std::pair<std::string, Header>(newKey,temp) );
	}
	else
	{
		mError(key_not_found);
	}
}

//Reset and validate all keys in headers and subheaders
void Document::revalidateAllKeys()
{
	//Reset head map keys
	this->resetKeys();
	
	//Loop through all headers
	for (auto &x: this->Head_Map)
	{
		x.second.resetKeys();
	}
}

//Adds a pair object to the map (with only strings)
void Document::addPair(std::string key, std::string val)
{
	this->Data_Map.addPair(key, val);
}

//Adds a pair object and asserts a type
void Document::addPair(std::string key, std::string val, int t)
{
	this->Data_Map.addPair(key, val, t);
}

//Set the name of the document
void Document::setName(std::string n)
{
	this->name = n;
}

//Set the alias of the document
void Document::setAlias(std::string a)
{
	this->alias = a;
}

//Set the name and alias of the document
void Document::setNameAliasPair(std::string n, std::string a, int s)
{
	this->SubHeader::setNameAliasPair(n, a, s);
}

//Set the state of the document
void Document::setState(int state)
{
	this->SubHeader::setState(state);
}

//Display the contents of the document
void Document::DisplayContents()
{
	std::cout << "Document Class\nName: " << this->getName();
	if (this->isAnchor())
		std::cout << " -> ANCHOR: " << this->getAlias() << "\n----------------------\n";
	else if (this->isAlias())
		std::cout << " -> ALIAS: " << this->getAlias() << "\n----------------------\n";
	else
		std::cout << "\n----------------------\n";
	std::cout << "Contents of this document include "<<this->size()<< " Headers and " <<this->Data_Map.size()<< " KeyValue Pairs...\n\n";
	
	//Display contents of KeyValueMap
	if (this->Data_Map.size() > 0)
		this->Data_Map.DisplayMap();
	
	//Loop through the SubHeaders
	for (auto &x: this->Head_Map)
	{
		x.second.DisplayContents();
	}
}

//Add a key to the Header map
void Document::addHeadKey(std::string key)
{
	Header temp(key);
	this->Head_Map.insert(this->Head_Map.end(), std::pair<std::string, Header> (key, temp) );
}

//Find and copy the contents of the Header Anchor to the Alias Header
void Document::copyAnchor2Alias(std::string alias, Header &ref)
{
	std::string key = ref.getName();
	ref = getAnchoredHeader(alias);
	ref.setNameAliasPair(key, alias, ALIAS);
	
}

//Return the size of the header map
int Document::size()
{
	return (int) this->Head_Map.size();
}

//Return the name of the document
std::string Document::getName()
{
	return this->name;
}

//Return the alias of the document
std::string Document::getAlias()
{
	return this->alias;
}

//Return the state of the document
int Document::getState()
{
	return this->state;
}

//Return true if document is Alias
bool Document::isAlias()
{
	if (this->getState() == ALIAS)
		return true;
	else
		return false;
}

//Return true if document is Anchor
bool Document::isAnchor()
{
	if (this->getState() == ANCHOR)
		return true;
	else
		return false;
}

//Return reference to the anchored header
Header& Document::getAnchoredHeader(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Head_Map)
	{
		if (x.second.isAnchor())
		{
			if (alias == x.second.getAlias())
			{
				foundit = true;
				current = x.first;
				break;
			}
		}
		else
		{
			current = x.first;
		}
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getHeader(current);
}

//Return reference to the header containing the Sub with the given alias
Header& Document::getHeadFromSubAlias(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Head_Map)
	{
		if (foundit == false)
		{
			for (auto &y: x.second.getSubMap())
			{
				if (y.second.isAnchor())
				{
					if (alias == y.second.getAlias())
					{
						foundit = true;
						current = x.first;
						break;
					}
				}
				else
				{
					current = x.first;
				}
			}
		}
		else
			break;
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getHeader(current);
}

/*
 =========================================== END Document Class ==============================================
 */

/*
 =========================================== YamlWrapper Class ==============================================
 */

//Default constructor
YamlWrapper::YamlWrapper()
{
}

//Default destructor
YamlWrapper::~YamlWrapper()
{
	Doc_Map.clear();
}

//Copy constructor
YamlWrapper::YamlWrapper(const YamlWrapper &yaml)
{
	Doc_Map = yaml.Doc_Map;
}

//Constructor by a document
YamlWrapper::YamlWrapper(std::string key, const Document &doc)
{
	Doc_Map.insert(Doc_Map.end(), std::pair<std::string, Document>(key, doc) );
}

//Equals operator overload
YamlWrapper& YamlWrapper::operator=(const YamlWrapper &yaml)
{
	if (this == &yaml)
		return *this;
	else
	{
		this->Doc_Map.clear();
		this->Doc_Map = yaml.Doc_Map;
		return *this;
	}
}

//Reference fetch overload
Document& YamlWrapper::operator()(const std::string key)
{
	try
	{
		return this->Doc_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Document>::iterator it;
		it = this->Doc_Map.end();
		return it->second;
	}
}

//Access operator overload
Document YamlWrapper::operator()(const std::string key) const
{
	try
	{
		return this->Doc_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Document>::const_iterator it;
		it = this->Doc_Map.end();
		return it->second;
	}
}

//Return the reference to the document map
std::map<std::string, Document> & YamlWrapper::getDocMap()
{
	return this->Doc_Map;
}

//Return reference to the document at the given key
Document& YamlWrapper::getDocument(std::string key)
{
	try
	{
		return this->Doc_Map.at(key);
	}
	catch (std::out_of_range)
	{
		throw;
		mError(key_not_found);
		std::map<std::string, Document>::iterator it;
		it = this->Doc_Map.end();
		return it->second;
	}
}

//Return a constant iterator for the end of the custom map
std::map<std::string, Document>::const_iterator YamlWrapper::end() const
{
	std::map<std::string, Document>::const_iterator it;
	it = this->Doc_Map.end();
	return it;
}

//Return an iterator for the end of the custom map
std::map<std::string, Document>::iterator YamlWrapper::end()
{
	std::map<std::string, Document>::iterator it;
	it = this->Doc_Map.end();
	return it;
}

//Return a constant iterator for the beginning of the custom map
std::map<std::string, Document>::const_iterator YamlWrapper::begin() const
{
	std::map<std::string, Document>::const_iterator it;
	it = this->Doc_Map.begin();
	return it;
}

//Return an iterator for the beginning of the custom map
std::map<std::string, Document>::iterator YamlWrapper::begin()
{
	std::map<std::string, Document>::iterator it;
	it = this->Doc_Map.begin();
	return it;
}

//Clear out memory allocated in YamlWrapper object
void YamlWrapper::clear()
{
	this->Doc_Map.clear();
}

//Reset all the keys for the document map
void YamlWrapper::resetKeys()
{
	std::string newKey;
	std::string oldKey;
	for (auto &x: this->Doc_Map)
	{
		oldKey = x.first;
		newKey = x.second.getName();
		std::swap(this->Doc_Map[newKey],x.second);
		this->Doc_Map.erase(x.first);
	}
}

//Change the given oldKey to the newKey in the map, if oldKey exists
void YamlWrapper::changeKey(std::string oldKey, std::string newKey)
{
	std::map<std::string, Document>::const_iterator it = this->Doc_Map.find(oldKey);
	if (it != this->Doc_Map.end())
	{
		Document temp(it->second);
		this->Doc_Map.erase(it);
		this->Doc_Map.insert(this->Doc_Map.end(), std::pair<std::string, Document>(newKey,temp) );
	}
	else
	{
		mError(key_not_found);
	}
}

//Revalidation of all keys and names in the structure
void YamlWrapper::revalidateAllKeys()
{
	//Reset keys in Document map
	this->resetKeys();
	
	//Loop through all document objects
	for (auto &x: this->Doc_Map)
	{
		//Reset keys in the Header maps contained within each Document Map
		x.second.revalidateAllKeys();
	}
}

//Display the contents of the YamlWrapper object
void YamlWrapper::DisplayContents()
{
	std::cout << "YamlWrapper Class\n";
	std::cout << "Contents of this YamlWrapper include "<<this->size()<< " Documents...\n\n";
	
	//Loop through the SubHeaders
	for (auto &x: this->Doc_Map)
	{
		std::cout << "-------------------------- START Document: " << x.second.getName() << " --------------------------------\n";
		x.second.DisplayContents();
		std::cout << "-------------------------- END Document: " << x.second.getName() << " --------------------------------\n\n";
	}
}

//Add a key to the document map
void YamlWrapper::addDocKey(std::string key)
{
	Document temp(key);
	this->Doc_Map.insert(this->Doc_Map.end(), std::pair<std::string, Document>(key, temp) );
}

//Find and copy the contents of the Document Anchor to the Alias Document
void YamlWrapper::copyAnchor2Alias(std::string alias, Document &ref)
{
	std::string key = ref.getName();
	ref = getAnchoredDoc(alias);
	ref.setNameAliasPair(key, alias, ALIAS);
	
}

//Return the size of the document map
int YamlWrapper::size()
{
	return (int) this->Doc_Map.size();
}

//Return reference to the Document in the map that is anchored with the given alias
Document& YamlWrapper::getAnchoredDoc(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Doc_Map)
	{
		if (x.second.isAnchor())
		{
			if (alias == x.second.getAlias())
			{
				foundit = true;
				current = x.first;
				break;
			}
		}
		else
		{
			current = x.first;
		}
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getDocument(current);
}

//Return reference to the document containing the header with the given alias
Document& YamlWrapper::getDocFromHeadAlias(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Doc_Map)
	{
		if (foundit == false)
		{
			for (auto &y: x.second.getHeadMap())
			{
				if (y.second.isAnchor())
				{
					if (alias == y.second.getAlias())
					{
						foundit = true;
						current = x.first;
						break;
					}
				}
				else
				{
					current = x.first;
				}
			}
		}
		else
			break;
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getDocument(current);
}

//Return reference to the document containing the subheader with the given alias
Document& YamlWrapper::getDocFromSubAlias(std::string alias)
{
	bool foundit = false;
	std::string current;
	for (auto &x: this->Doc_Map)
	{
		if (foundit == false)
		{
			for (auto &y: x.second.getHeadMap())
			{
				if (foundit == false)
				{
					for (auto &z: y.second.getSubMap())
					{
						if (z.second.isAnchor())
						{
							if (alias == z.second.getAlias())
							{
								foundit = true;
								current = x.first;
								break;
							}
						}
						else
						{
							current = x.first;
						}
					}
				}
				else
					break;
			}
		}
		else
			break;
	}
	if (foundit == false)
		mError(anchor_alias_dne);
	return this->getDocument(current);
}

/*
 =========================================== END YamlWrapper Class ==============================================
 */

/*
 =========================================== START yaml_cpp_class ==============================================
 */

//Default constructor
yaml_cpp_class::yaml_cpp_class()
{
}

//Default destructor
yaml_cpp_class::~yaml_cpp_class()
{
	yaml_wrapper.clear();
}

//Set the input file for reading
int yaml_cpp_class::setInputFile(const char *file)
{
	int success = 0;
	
	this->file_name = file;
	if (this->file_name == NULL || *this->file_name == ' ')
	{
		mError(file_dne);
		return -1;
	}
	this->input_file = fopen(this->file_name, "r");
	if(!yaml_parser_initialize(&this->token_parser))
	{
		mError(initial_error);
		return -1;
	}
	if (this->input_file == NULL)
	{
		mError(file_dne);
		return -1;
	}
	yaml_parser_set_input_file(&this->token_parser, this->input_file);
	
	return success;
}

//Reads the input file and store the results into the YamlWrapper
int yaml_cpp_class::readInputFile()
{
	int success = 0;
	
	bool new_block = false;
	bool block_key = false;
	bool sub_block = false;
	bool doc_start = false;
	std::string current_key;
	std::string current_value;
	std::string current_head;
	std::string current_sub;
	std::string current_doc;
	std::string current_alias;
	std::string current_anchor;
	
	do
	{
		yaml_parser_scan(&this->token_parser, &this->current_token);
		
		//Switch case that is used to change code behavior based on token type
		switch(this->current_token.type)
		{
			case YAML_STREAM_START_TOKEN:
			{
				break;
			}
			case YAML_STREAM_END_TOKEN:
			{
				break;
			}
			case YAML_KEY_TOKEN:
			{
				//Check to see if we are in a header/block
				if (this->previous_token.type == YAML_BLOCK_MAPPING_START_TOKEN)
				{
					block_key = true;
				}
				else
				{
					block_key = false;
				}
				break;
			}
			case YAML_VALUE_TOKEN:
			{
				break;
			}
			case YAML_BLOCK_SEQUENCE_START_TOKEN:
			{
				break;
			}
			case YAML_BLOCK_ENTRY_TOKEN:
			{
				//Check to see if we are in a subheader or main header
				if (new_block == true)
					sub_block = true;
				else
					sub_block = false;
				new_block = true;
				break;
			}
			case YAML_BLOCK_END_TOKEN:
			{
				if (sub_block == true)
					sub_block = false;
				else
					new_block = false;
				break;
			}
			case YAML_BLOCK_MAPPING_START_TOKEN:
			{
				break;
			}
			case YAML_SCALAR_TOKEN:
			{
				//Check to see what the type is for the current scalar
				//NOTE: Current scalar is always the type of the previous token
				if (this->previous_token.type == YAML_KEY_TOKEN)
				{
					current_key = (char*) this->current_token.data.scalar.value;
					
					//Check to see what kind of key needs to be added
					if (doc_start == false)
					{
						current_doc = current_key;
						this->yaml_wrapper.addDocKey(current_doc);
						this->yaml_wrapper(current_doc).setName(current_doc);
					}
					else if (block_key == true && sub_block == false && new_block == true)
					{
						current_head = current_key;
						this->yaml_wrapper(current_doc).addHeadKey(current_head);
						this->yaml_wrapper(current_doc)(current_head).setName(current_head);
					}
					else if (block_key == true && sub_block == true &&  new_block == true)
					{
						current_sub = current_key;
						this->yaml_wrapper(current_doc)(current_head).addSubKey(current_sub);
						this->yaml_wrapper(current_doc)(current_head)(current_sub).setName(current_sub);
					}
					else
					{
						//Check to see if we are in a document, header, or sub-header
						if (sub_block == false && new_block == false)
						{
							this->yaml_wrapper(current_doc).getDataMap().addKey(current_key);
						}
						else if (block_key == false && sub_block == false && new_block == true)
						{
							this->yaml_wrapper(current_doc)(current_head).getDataMap().addKey(current_key);
						}
						else if (block_key == false && sub_block == true && new_block == true)
						{
							this->yaml_wrapper(current_doc)(current_head)(current_sub).getMap().addKey(current_key);
						}
						else
						{
							mError(not_a_token);
							return -1;
						}
					}
					
				}
				else if (this->previous_token.type == YAML_VALUE_TOKEN && doc_start == true)
				{
					current_value = (char*) this->current_token.data.scalar.value;
					//Check to see if we are in a document, header, or sub-header
					if (sub_block == false && new_block == false)
					{
						this->yaml_wrapper(current_doc).getDataMap().editValue4Key(current_value, current_key);
					}
					else if (block_key == false && sub_block == false && new_block == true)
					{
						this->yaml_wrapper(current_doc)(current_head).getDataMap().editValue4Key(current_value, current_key);
					}
					else if (block_key == false && sub_block == true && new_block == true)
					{
						this->yaml_wrapper(current_doc)(current_head)(current_sub).getMap().editValue4Key(current_value, current_key);
					}
					else
					{
						mError(not_a_token);
						return -1;
					}
	
				}
				else
				{
					mError(not_a_token);
					return -1;
				}
				break;
			}
			case YAML_NO_TOKEN:
			{
				mError(not_a_token);
				return -1;
			}
			case YAML_VERSION_DIRECTIVE_TOKEN:
			{
				break;
			}
			case YAML_TAG_DIRECTIVE_TOKEN:
			{
				break;
			}
			case YAML_DOCUMENT_START_TOKEN:
			{
				doc_start = true;
				break;
			}
			case YAML_DOCUMENT_END_TOKEN:
			{
				doc_start = false;
				break;
			}
			case YAML_FLOW_SEQUENCE_START_TOKEN:
			{
				break;
			}
			case YAML_FLOW_SEQUENCE_END_TOKEN:
			{
				break;
			}
			case YAML_FLOW_MAPPING_START_TOKEN:
			{
				break;
			}
			case YAML_FLOW_MAPPING_END_TOKEN:
			{
				break;
			}
			case YAML_FLOW_ENTRY_TOKEN:
			{
				break;
			}
			case YAML_ALIAS_TOKEN:
			{
				current_alias = (char *) this->current_token.data.alias.value;
				if (doc_start == true && sub_block == false && new_block == false)
				{
					mError(not_a_token);
					return -1;
				}
				
				//Check to see what it is an alias of
				if (doc_start == false)
				{
					// Search for document anchor and copy contents
					this->yaml_wrapper.copyAnchor2Alias(current_alias, this->yaml_wrapper(current_doc));
				}
				else if (block_key == true && sub_block == false && new_block == true)
				{
					// Search through all documents for header anchor and copy contents
					this->yaml_wrapper.getDocFromHeadAlias(current_alias).copyAnchor2Alias(current_alias, this->yaml_wrapper(current_doc)(current_head));
				}
				else if (block_key == true && sub_block == true &&  new_block == true)
				{
					// Search through all documents and headers for sub-header anchor and copy contents
					this->yaml_wrapper.getDocFromSubAlias(current_alias).getHeadFromSubAlias(current_alias).copyAnchor2Alias(current_alias, this->yaml_wrapper(current_doc)(current_head)(current_sub));
				}
				
				break;
			}
			case YAML_ANCHOR_TOKEN:
			{
				current_anchor = (char *) this->current_token.data.anchor.value;
				if (doc_start == true && sub_block == false && new_block == false)
				{
					mError(not_a_token);
					return -1;
				}
				
				//Check to see what it is an anchor of
				if (doc_start == false)
				{
					this->yaml_wrapper(current_doc).setAlias(current_anchor);
					this->yaml_wrapper(current_doc).setState(ANCHOR);
				}
				else if (block_key == true && sub_block == false && new_block == true)
				{
					this->yaml_wrapper(current_doc)(current_head).setAlias(current_anchor);
					this->yaml_wrapper(current_doc)(current_head).setState(ANCHOR);
				}
				else if (block_key == true && sub_block == true &&  new_block == true)
				{
					this->yaml_wrapper(current_doc)(current_head)(current_sub).setAlias(current_anchor);
					this->yaml_wrapper(current_doc)(current_head)(current_sub).setState(ANCHOR);
				}
				
				break;
			}
			case YAML_TAG_TOKEN:
			{
				break;
			}
			default:
			{
				mError(not_a_token);
				return -1;
			}
				
		}//END SWITCH CASE
		this->previous_token = this->current_token;
		
	} while (this->current_token.type != YAML_STREAM_END_TOKEN);
	
	return success;
}

//Deletes the yaml_c objects and closes the input file
int yaml_cpp_class::cleanup()
{
	int success = 0;
	
	yaml_token_delete(&this->previous_token);
	yaml_token_delete(&this->current_token);
	yaml_parser_delete(&this->token_parser);
	fclose(this->input_file);
	
	return success;
}

//Runs full execution of object intialization, data reading, and cleaning of the structure
int yaml_cpp_class::executeYamlRead(const char *file)
{
	int success = 0;
	
	//Try to open file and initialize yaml objects
	success = this->setInputFile(file);
	if (success != 0)
	{
		mError(read_error);
		return -1;
	}
	
	//Try to read the opened file and store the data
	success = this->readInputFile();
	if (success != 0)
	{
		mError(read_error);
		return -1;
	}
	
	//Close the file and delete temporary objects
	success = this->cleanup();
	
	return success;
}

//Display contents of the read to the console window
void yaml_cpp_class::DisplayContents()
{
	this->yaml_wrapper.DisplayContents();
}

//Return reference to the YamlWrapper object
YamlWrapper& yaml_cpp_class::getYamlWrapper()
{
	return this->yaml_wrapper;
}

/*
 =========================================== END yaml_cpp_class ==============================================
 */

// Convert input to all lower case
std::string allLower(const std::string &input)
{
	std::string copy = input;
	for (int i=0; i<copy.size(); i++)
		copy[i] = tolower(copy[i]);
	return copy;
}

//Return true if arg is even
bool isEven(int n)
{
	if (n%2 == 0)
		return true;
	else
		return false;
}

int YAML_WRAPPER_TESTS()
{
	int success = 0;
	
	ValueTypePair test1("apple",STRING);
	ValueTypePair test2;
	
	test2.editPair("-1e12", DOUBLE);
	
	test1.DisplayPair();
	test2.DisplayPair();
	
	std::string num = "-1.56E-23";
	double val = atof(num.c_str());
	std::cout << num << " = " << val << std::endl;
	
	test2.findType();
	
	ValueTypePair test3("3",INT);
	test3.findType();
	test3.DisplayPair();
	test3.assertType(STRING);
	test3.DisplayPair();
	
	KeyValueMap maptest1("key1","value1");
	
	ValueTypePair val1;
	
	val1 = maptest1["key1"];
	val1.DisplayPair();
	maptest1.DisplayMap();
	maptest1["key1"] = test3;
	maptest1.addPair("key2", test1);
	maptest1.addPair("key3", test2);
	maptest1.addPair("key4", "value4");
	maptest1.addPair("key5", "2.34E-12");
	maptest1.addPair("key6", "FALSE");
	maptest1.addPair("key7", "132");
	maptest1.addPair("key8", "35.4",INT);
	maptest1.DisplayMap();
	
	
	val1.DisplayPair();
	maptest1.findType("key1");
	maptest1.assertType("key8", STRING);
	
	maptest1.findAllTypes();
	maptest1.DisplayMap();
	
	ValueTypePair test4 = maptest1.getPair("key5");
	test4.DisplayPair();
	
	maptest1.addKey("key9");
	maptest1.editValue4Key("1.0", "key9");
	maptest1.editValue4Key("1.0", INT, "key9");
	maptest1.DisplayMap();
	
	SubHeader sub1("Sub1",maptest1);
	
	sub1["key8"].DisplayPair();
	test4 = sub1["key1"];
	test4.DisplayPair();
	
	sub1.getMap().DisplayMap();
	
	std::cout << "key5 (" << sub1.getMap().getType("key5") << ") = " << sub1.getMap().getDouble("key5") << std::endl;
	std::cout << "key9 (" << sub1.getMap().getType("key9") << ") = " << sub1.getMap().getInt("key9") << std::endl;
	
	sub1.setNameAliasPair("apples", "fruit", ALIAS);
	sub1.DisplayContents();
	
	SubHeader sub2;
	
	sub2.setNameAliasPair("ShoppingList", " " , NONE);
	
	sub2.getMap().addKey("banana");
	sub2.getMap().addKey("gravy");
	sub2.getMap().addKey("screw");
	sub2.getMap().addKey("cherry");
	
	sub2.getMap().editValue4Key("12", "banana");
	sub2.getMap().editValue4Key("10 oz", "gravy");
	sub2.getMap().editValue4Key("six", "screw");
	sub2.getMap().editValue4Key("22.5 lbs", "cherry");
	
	sub2.getMap().findAllTypes();
	
	sub2.DisplayContents();
	
	std::cout << "lbs of cherries = " << sub2.getMap().getDouble("cherry") << std::endl;
	
	
	Header head1("ShopList",sub2.getMap());
	
	head1["banana"].DisplayPair();
	test4 = head1["cherry"];
	test4.DisplayPair();
	
	Header head2("keylist",sub1);
	
	head2("keylist").DisplayContents();
	SubHeader newsub;
	newsub = head2("keylist");
	newsub.DisplayContents();
	
	head2.setNameAliasPair("the keylist", "none", NONE);
	head2.getDataMap().addPair("Monday", "29");
	head2.getDataMap().addPair("Tuesday", "30");
	head2.getDataMap().addPair("Wednesday", "31");
	
	head2.addSubKey("NextSub");
	head2("NextSub") = sub2;
	
	//head2.resetKeys();
	head2.changeKey("NextSub", "ShoppingList");
	
	//head2("ShoppingList").DisplayContents();
	head2.getSubHeader("ShoppingList").getMap().addPair("apples", "27");
	head2.getSubHeader("ShoppingList").getMap().addPair("pie", "2.7 inch");
	
	head2.addSubKey("NextList");
	head2.getSubHeader("NextList").getMap().addPair("test1", "this is a test");
	head2.getSubHeader("NextList").getMap().addPair("test2", "FALSE");
	head2.getSubHeader("NextList").getMap().addPair("test3", "1.4e-12");
	head2.getSubHeader("NextList").getMap().addPair("test4", "23");
	
	head2.setAlias("head2");
	head2.setState(ANCHOR);
	
	head2.DisplayContents();
	
	//Build out document
	Document doc1;
	
	doc1.setNameAliasPair("Document1", "DocTest", ALIAS);
	doc1.getDataMap().addPair("Key1", "26.5");
	doc1.getDataMap().addPair("Key2", "3");
	doc1.getDataMap().addPair("Key3", "true");
	doc1.getDataMap().addPair("Key4", "string test");
	
	doc1.addHeadKey("Header1");
	doc1.getHeader("Header1").getDataMap().addPair("HeadKey1", "header key test");
	doc1.getHeader("Header1").getDataMap().addPair("Deadkey2", "dead key?");
	
	doc1.getHeader("Header1").addSubKey("SubHead1");
	doc1.getHeader("Header1").getSubHeader("SubHead1").addPair("maples", "syrup");
	doc1.getHeader("Header1").getSubHeader("SubHead1").addPair("choco", "latte");
	doc1.getHeader("Header1").addSubKey("SubHead2");
	doc1.getHeader("Header1").getSubHeader("SubHead2").getMap().addPair("thursday?", "false");
	
	doc1.addHeadKey("Header2");
	doc1.getHeader("Header2").getDataMap().addPair("how deep", "does the rabbit hole go?");
	doc1.getHeader("Header2").addPair("about", "36 ft");
	
	doc1.getHeader("Header2").addSubKey("SubHead1");
	doc1.getHeader("Header2").getSubHeader("SubHead1").getMap().addPair("that's right", "this sub is unique");
	doc1.getHeader("Header2").addSubKey("SubHead2");
	doc1.getHeader("Header2").getSubHeader("SubHead2").addPair("isUnique", "true");
	
	
	doc1.DisplayContents();
	
	std::cout << "how many feet? " << doc1.getHeader("Header2").getDataMap().getInt("about") << std::endl << std::endl;
	
	std::cout << "--------------- FINAL TEST before Wrapping the C Parser -----------------\n\n";
	
	YamlWrapper yaml;
	
	yaml.addDocKey("Doc1");
	yaml.getDocument("Doc1").setAlias("TestAnchor");
	yaml.getDocument("Doc1").setState(ANCHOR);
	yaml.getDocument("Doc1").addPair("doc1-key1", "273.15");
	yaml.getDocument("Doc1").addPair("doc1-key2", "3");
	yaml.getDocument("Doc1").addHeadKey("Doc1-Head1");
	yaml.getDocument("Doc1").getHeader("Doc1-Head1").addPair("doc1-head1-key1", "true");
	yaml.getDocument("Doc1").addHeadKey("Doc1-Head2");
	yaml.getDocument("Doc1").getHeader("Doc1-Head2").addSubKey("Doc1-Head2-Sub1");
	yaml.getDocument("Doc1").getHeader("Doc1-Head2").getSubHeader("Doc1-Head2-Sub1").setAlias("that sub", ANCHOR);
	yaml.getDocument("Doc1").getHeader("Doc1-Head2").getSubHeader("Doc1-Head2-Sub1").addPair("doc1-head2-sub1-key1", "lots of space");
	
	yaml.addDocKey("Doc2");
	yaml.getDocument("Doc2").setAlias("TestAnchor");
	yaml.getDocument("Doc2").setState(ALIAS);
	yaml.copyAnchor2Alias(yaml.getDocument("Doc2").getAlias(), yaml.getDocument("Doc2"));
	//yaml.getAnchoredDoc("TestAnchor").DisplayContents();
	yaml.getDocument("Doc2").addHeadKey("Doc2-Head1");
	yaml.getDocument("Doc2").getHeader("Doc2-Head1").setAlias("HeadAnchor");
	yaml.getDocument("Doc2").getHeader("Doc2-Head1").setState(ANCHOR);
	yaml.getDocument("Doc2").getHeader("Doc2-Head1").addPair("added stuff", "true");
	//Note: we can add stuff to the aliased object after copying contents from the anchor
	yaml.getDocument("Doc2").getHeader("Doc1-Head1").getDataMap().editValue4Key("false", "doc1-head1-key1");
	//Note: we can also change the existing data copied over from the anchored-aliased pair
	//		The edited value is only changed in the alias, the anchor remains constant.
	
	
	yaml.addDocKey("Doc3");
	yaml.getDocument("Doc3").addHeadKey("Header-Ref");
	yaml.getDocument("Doc3").getHeader("Header-Ref").setAlias("HeadAnchor");
	yaml.getDocument("Doc3").getHeader("Header-Ref").setState(ALIAS);
	
	//Copy from Doc2 to Doc3; finding the Alias/Anchor "HeadAnchor" and copying into "Header-Ref" of Doc3
	//yaml.getDocument("Doc2").copyAnchor2Alias(yaml.getDocument("Doc3").getHeader("Header-Ref").getAlias(), yaml.getDocument("Doc3").getHeader("Header-Ref"));
	//Problem:: I had to know ahead of time what document the anchor header was in
	
	yaml.getDocFromHeadAlias(yaml.getDocument("Doc3").getHeader("Header-Ref").getAlias()).copyAnchor2Alias(yaml.getDocument("Doc3").getHeader("Header-Ref").getAlias(), yaml.getDocument("Doc3").getHeader("Header-Ref"));
	//This function will perform the copy implicitly, but requires alot of searching
	//Function only works if anchor has been previously established
	//There cannot be two anchors with the same name, but there can be many aliases of the same name
	
	yaml.getDocument("Doc3").getHeader("Header-Ref").addSubKey("Sub-Ref");
	yaml.getDocument("Doc3").getHeader("Header-Ref").getSubHeader("Sub-Ref").setAlias("that sub", ALIAS);
	
	//This entire line of code is used to find the subheader anchor and copy it's contents into the alias specified from above
	yaml.getDocFromSubAlias(yaml.getDocument("Doc3").getHeader("Header-Ref").getSubHeader("Sub-Ref").getAlias()).getHeadFromSubAlias(yaml.getDocument("Doc3").getHeader("Header-Ref").getSubHeader("Sub-Ref").getAlias()).copyAnchor2Alias(yaml.getDocument("Doc3").getHeader("Header-Ref").getSubHeader("Sub-Ref").getAlias(), yaml.getDocument("Doc3").getHeader("Header-Ref").getSubHeader("Sub-Ref"));
	
	//All functions are working as they should, the YamlWrapper object is complete to this point
	//Now move on to the read parser and store the c parser read into the object.
	
	yaml.DisplayContents();
	
	std::cout << "--------------------------- END FINAL TEST ------------------------------\n";
	
	
	return success;
}

//Run tests for YAML_CPP_CLASS
int YAML_CPP_TEST(const char *file)
{
	int success = 0;

	yaml_cpp_class yaml_object;
	success = yaml_object.executeYamlRead(file);
	if (success != 0)
	{
		mError(simulation_fail);
		return -1;
	}
	yaml_object.DisplayContents();
	
	return success;
}
