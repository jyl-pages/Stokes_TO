//////////////////////////////////////////////////////////////////////////
// Hashtable
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Hashtable_h__
#define __Hashtable_h__
#include <unordered_map>
#include <unordered_set>
#include "Common.h"

////Hash key instantiations
namespace std{
template<> struct hash<Vector2i>
{typedef Vector2i argument_type;typedef std::size_t result_type;
	result_type operator()(argument_type const& arg) const
	{result_type const h1(std::hash<int>()(arg[0]));result_type const h2(std::hash<int>()(arg[1]));return h1^(h2<<1);}
};
template<> struct hash<Vector3i>
{typedef Vector3i argument_type;typedef std::size_t result_type;
	result_type operator()(argument_type const& arg) const
	{result_type const h1(std::hash<int>()(arg[0]));result_type const h2(std::hash<int>()(arg[1]));
	result_type const h3(std::hash<int>()(arg[2]));return h1^(h2<<1)^h3;}
};
template<> struct hash<Vector4i>
{typedef Vector4i argument_type;typedef std::size_t result_type;
	result_type operator()(argument_type const& arg) const
	{result_type const h1(std::hash<int>()(arg[0]));result_type const h2(std::hash<int>()(arg[1]));
	result_type const h3(std::hash<int>()(arg[2]));result_type const h4(std::hash<int>()(arg[3]));return h1^(h2<<1)^h3^(h4<<2);}
};}

////Hashtable
template<class T_KEY,class T> using Hashtable=std::unordered_map<T_KEY,T>;
template<class T_KEY,class T> using HashtableMultiValue=std::unordered_multimap<T_KEY,T>;
////Hashset
template<class T_KEY> using Hashset=std::unordered_set<T_KEY>;

////Hashtable operations
template<class T_KEY,class T> bool Has(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const T& v)
{auto range=hashtable.equal_range(key);
for(auto iter=range.first;iter!=range.second;iter++){if(iter->second==v){return true;}}return false;}

template<class T_KEY> bool Has(Hashset<T_KEY>& hashset,const T_KEY& key)
{auto iter=hashset.find(key);return iter!=hashset.end();}

template<class T_KEY,class T> bool Add(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const T& v,bool check_duplicate=false)
{bool add=!check_duplicate||!Has(hashtable,key,v);if(add)hashtable.insert(std::make_pair(key,v));return add;}

template<class T_KEY,class T> bool Add(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const Array<T>& values,bool check_duplicate=false)
{bool add=true;for(auto& v:values)if(!Add(hashtable,key,v,check_duplicate))add=false;return add;}

template<class T_KEY> bool Add(Hashset<T_KEY>& hashset,const T_KEY& key)
{auto iter=hashset.insert(key);return iter.second;}

template<class T_KEY,class T> bool Remove(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const T& v)
{auto range=hashtable.equal_range(key);for(auto iter=range.first;iter!=range.second;iter++){if(iter->second==v){hashtable.erase(iter);return true;}}return false;}

template<class T_KEY,class T> int Remove(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key)
{return (int)(hashtable.erase(key));}

template<class T_KEY> bool Remove(Hashset<T_KEY>& hashset,const T_KEY& key)
{int n=(int)hashset.erase(key);return n!=0;}

template<class T_KEY,class T> bool Remove_Multiple(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const Array<T>& values)
{for(const auto& v:values)return Remove(hashtable,key,v);}

template<class T_KEY,class T> bool Replace(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const T& v0,const T& v1)
{if(v0==v1)return true;auto range=hashtable.equal_range(key);
for(auto iter=range.first;iter!=range.second;iter++){if(iter->second==v0){iter->second=v1;return true;}}return false;}

template<class T_KEY,class T> bool Replace_Or_Add(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,const T& v0,const T& v1)
{auto range=hashtable.equal_range(key);bool replaced=false;
for(auto iter=range.first;iter!=range.second;iter++){if(iter->second==v0){iter->second=v1;replaced=true;break;}}
if(!replaced)Add(hashtable,key,v1);return true;}

template<class T_KEY,class T> void Key_List(const HashtableMultiValue<T_KEY,T>& hashtable,List<T_KEY>& keys)
{for(const auto& iter:hashtable){keys.push_back(iter.first);}}

template<class T_KEY,class T> void Value_List(const HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,List<T>& values)
{auto range=hashtable.equal_range(key);for(auto iter=range.first;iter!=range.second;iter++){values.push_back(iter->second);}}
template<class T_KEY,class T> void Value_Array(const HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key,Array<T>& values)
{auto range=hashtable.equal_range(key);for(auto iter=range.first;iter!=range.second;iter++){values.push_back(iter->second);}}
template<class T_KEY,class T> int Value_Count(const HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key)
{return (int)hashtable.count(key);}

template<class T_KEY,class T> void Replace_Key(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key_old,const T_KEY& key_new)
{List<T> values;Value_List(hashtable,key_old,values);Remove(hashtable,key_old);for(const auto& v:values)Add(hashtable,key_new,v);}

template<class T_KEY,class T> void Merge_Keys(HashtableMultiValue<T_KEY,T>& hashtable,const T_KEY& key_old_0,const T_KEY& key_old_1,const T_KEY& key_new)
{List<T> v0;Value_List(hashtable,key_old_0,v0);List<T> v1;Value_List(hashtable,key_old_1,v1);
for(const auto& v:v1){if(std::find(v0.begin(),v0.end(),v)==v0.end())v0.push_back(v);}
Remove(hashtable,key_old_0);Remove(hashtable,key_old_1);for(const auto& v:v0)Add(hashtable,key_new,v);}

template<class T,class T_FUNC> void Heap_To_Hashset(const Heap<T,T_FUNC>& heap,Hashset<T>& hashset)
{Heap<T,T_FUNC> copy=heap;while(!copy.empty()){hashset.insert(copy.top());copy.pop();}}

template<class T> void Set_Difference(List<T>& a,List<T>& b,List<T>& dif,const bool sorted=false)
{if(!sorted){a.sort();b.sort();}std::set_difference(a.begin(),a.end(),b.begin(),b.end(),std::inserter(dif,dif.begin()));}

//////////////////////////////////////////////////////////////////////////
////sorting related functions

////Vector index sorting
inline Vector1i Sorted(const Vector1i& v){return v;}
inline Vector2i Sorted(const Vector2i& v){return v[0]>v[1]?v:Vector2i(v[1],v[0]);}
inline Vector3i Sorted(const Vector3i& _v){Vector3i v=_v;if(v[0]>v[1])std::swap(v[0],v[1]);if(v[0]>v[2])std::swap(v[0],v[2]);if(v[1]>v[2])std::swap(v[1],v[2]);return v;}
////Vector index sorting by preserving cw/ccw order
inline Vector<int,2> Reordered(const Vector<int,2>& e,const int v0)
{return v0==e[0]?Vector<int,2>(e[0],e[1]):Vector<int,2>(e[1],e[0]);}
inline Vector<int,3> Reordered(const Vector<int,3>& e,const int v0)
{return v0==e[0]?Vector<int,3>(e[0],e[1],e[2]):v0==e[1]?Vector<int,3>(e[1],e[2],e[0]):Vector<int,3>(e[2],e[0],e[1]);}
inline Vector<int,4> Reordered(const Vector<int,4>& e,const int v0)
{return v0==e[0]?Vector<int,4>(e[0],e[1],e[2],e[3]):v0==e[1]?Vector<int,4>(e[1],e[2],e[0],e[3]):v0==e[2]?Vector<int,4>(e[2],e[0],e[1],e[3]):Vector<int,4>(e[3],e[1],e[0],e[2]);}

////Vector index sorting for a unique key value
inline Vector<int,2> Unique_Ordered(const Vector<int,2>& e){return Sorted(e);}
inline Vector<int,3> Unique_Ordered(const Vector<int,3>& e){return Reordered(e,e.minCoeff());}
inline Vector<int,4> Unique_Ordered(const Vector<int,4>& e)
{Vector<int,4> v1=Reordered(e,e.minCoeff());Vector<int,3> tri=Unique_Ordered(Vector<int,3>(v1[1],v1[2],v1[3]));v1[1]=tri[0];v1[2]=tri[1];v1[3]=tri[2];return v1;}
inline Vector<int,2> Unique_Ordered(const int v0,const int v1){return Unique_Ordered(Vector2i(v0,v1));}
inline Vector<int,3> Unique_Ordered(const int v0,const int v1,const int v2){return Unique_Ordered(Vector3i(v0,v1,v2));}
inline Vector<int,4> Unique_Ordered(const int v0,const int v1,const int v2,const int v3){return Unique_Ordered(Vector4i(v0,v1,v2,v3));}

////Vector index sorting for a reversely-ordered key value
inline Vector<int,2> Reversed(const Vector<int,2>& e){return Vector<int,2>(e[1],e[0]);}
inline Vector<int,3> Reversed(const Vector<int,3>& e){return Vector<int,3>(e[2],e[1],e[0]);}
inline Vector<int,4> Reversed(const Vector<int,4>& e){return Vector<int,4>(e[0],e[1],e[3],e[2]);}
template<int d> inline Vector<int,d> Unique_Reversed(const Vector<int,d>& e){return Unique_Ordered(Reversed(e));}	////2D is not reversed

////Vector index key comparison
inline bool Is_Same_Order(const Vector<int,3>& e1,const Vector<int,3>& e2){return e1==Reordered(e2,e1[0]);}
inline bool Is_Same_Order(const Vector<int,2>& edge,const Vector<int,3>& tri){return Reordered(tri,edge[0])[1]==edge[1];}
#endif