//////////////////////////////////////////////////////////////////////////
// Timer cpp
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Timer.h"

void TimerRecord::Add(double cur_time)
{
	total += cur_time;
	cur = cur_time;
	num++;
}

std::pair<double, double> TimerRecord::Profile(void)
{
	return std::pair<double, double>(cur, total / num);
}


template<class T>
void Timer<T>::Reset()
{
	start=std::chrono::system_clock::now();
}

template<class T>
T Timer<T>::Elapse()
{
	std::chrono::time_point<std::chrono::system_clock> end=std::chrono::system_clock::now();
	std::chrono::duration<T,std::ratio<1,1000> > elapse=end-start;
	return elapse.count();
}

template<class T>
T Timer<T>::Elapse_And_Reset()
{
	T elapse=Elapse();
	Reset();
	return elapse;
}

template<class T>
void Timer<T>::Elapse_And_Output(const std::string& message)
{
	T elapse=Elapse();
	std::cout << "Timer for " << message << ": " << elapse << " ms" << std::endl;
}

template<class T>
void Timer<T>::Elapse_And_Output_And_Reset(const std::string& message)
{
	Elapse_And_Output(message);
	Reset();
}


template<class T>
inline void Timer<T>::Begin_Loop(void)
{
	Reset();
}

template<class T>
void Timer<T>::Record(const std::string& name)
{
	int k = -1;
	auto iter = name_index.find(name);
	if (iter == name_index.end()) {
		k = records.size();
		name_index[name] = k;
		records.push_back(TimerRecord(name));//push a zero record
	}
	else {
		k = iter->second;
	}
	records[k].Add(Elapse_And_Reset() / 1000.0);//Elapse_And_Reset() returns in ms
}

template<class T>
void Timer<T>::End_Loop_And_Output(std::ostream& out)
{
	out << "Time record cur(avg) in seconds: ";
	for (int i = 0; i < records.size(); i++) {
		std::pair<double, double> rec = records[i].Profile();
		out << records[i].name << ": " << rec.first << "(" << rec.second << "), ";
	}
	out << "\n";
}

template class Timer<double>;