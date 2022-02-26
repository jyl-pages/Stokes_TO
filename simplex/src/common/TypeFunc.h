//////////////////////////////////////////////////////////////////////////
// Type function
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __AuxTypeFunc_h__
#define __AuxTypeFunc_h__

#define DeclareCheckFunc(class_name,func_name)\
template<class TYPE> class CheckFunc {\
private:template<class T> static decltype(std::declval<T>().f(),std::true_type()) check(int);\
template<class T> static std::false_type check(...);\
public:static constexpr bool value=std::is_same<decltype(check<real>(0)),std::true_type>::value;};

template<bool condition,class T1,class T2> struct If{typedef T1 Type;};
template<class T1,class T2> struct If<false,T1,T2>{typedef T2 Type;};

#endif