/*
 * Luis Cruz
 * Center for Polymer Studies
 * Boston University
 * September, 1996
 * ccruz@bu.edu
 *
 * File: 	infobjarray.h
 * Purpose:
 *    Template for a list of objects addresable by index number.
 *    The objects must:
 *	1. know how to initialize themselves when 'newed'
 * 	2. how to print themselves.
 * 	3. have an assignment operator.
 * 	4. can compare two using '=='.
 * 	5. have a 'clear' that puts to zero.
 *
 */

#ifndef INFOBJARRAY_H
#define INFOBJARRAY_H

#include <iostream>

#include <stdlib.h>

using namespace std;

template <class T, int MAX>
class infobjarray
{
	friend ostream &operator<<(ostream &os, infobjarray<T,MAX> &z)
	{
		infobjarray<T,MAX>* ptr_cout = &z;

		int i = 0 ;
		int j = 0 ;
		while ( i < z.tot_dim )
		{
			os << ptr_cout->data[j] << "\n";
			i++;
			j++;
			if ( j >= MAX )
			{
			    if ( ptr_cout->ptr_next != NULL )
			    {
				ptr_cout = ptr_cout->ptr_next;
				j -= MAX;
			    }
			}
		}
		return os;
	}        

public:
	infobjarray(void)
	{
		ptr_next = NULL;
		tot_dim = MAX;
		last_one = -1;
		data = new T[MAX];
	}       
	~infobjarray(void)
	{
		delete ptr_next;
		ptr_next = NULL;

		delete [] data;
		data = NULL;
	}       
private:
	infobjarray(infobjarray<T,MAX>& v);
	infobjarray &operator=(infobjarray<T,MAX>& v);
public:
	inline void add_one(const T& one)	
		// ADD AN ELEMENT OF CLASS T AT THE END.
	{
		last_one++;
		(*this)[last_one] = one;
	}       
	inline int erase_one(T& one)
		// RETURN 1, FOUND IT.
		// RETURN 0, NOT FOUND IT.
	{
		infobjarray<T,MAX>* ptr_erase = this;

		int i = 0 ;
		int j = 0 ;
		while ( i < this->last_one )
		{        
			if ( ptr_erase->data[j] == one )
			{
				ptr_erase->data[j] = (*this)[last_one];
				(*this)[last_one].clear();
				last_one--;
				return 1;
			}
			else
			{
				i++;
				j++;
				if ( j >= MAX )
				{
				     if ( ptr_erase->ptr_next != NULL )
				     {
					 ptr_erase = ptr_erase->ptr_next;
					 j -= MAX;
				     }    
				}
			}
		}        
		return 0;
	}       
	inline int erase_one_and_collapse(T& one)
		// RETURN 1, FOUND IT.
		// RETURN 0, NOT FOUND IT.
	{
		infobjarray<T,MAX>* ptr_erase = this;

		int i = 0 ;
		int j = 0 ;
		while ( i < this->last_one )
		{        
			if ( ptr_erase->data[j] == one )
			{
				for (int run=i; run<last_one; run++)
				{
				    (*this)[run] = (*this)[run+1];
				}
				(*this)[last_one].clear();
				last_one--;
				return 1;
			}
			else
			{
				i++;
				j++;
				if ( j >= MAX )
				{
				     if ( ptr_erase->ptr_next != NULL )
				     {
					 ptr_erase = ptr_erase->ptr_next;
					 j -= MAX;
				     }    
				}
			}
		}        
		return 0;
	}       
	inline void erase_last_one(void)
	{
		if (last_one < 0) return;

		(*this)[last_one].clear();
		last_one--;
	}       
	inline int is_not_empty(void)
	{
		if (last_one<0)
			return 0;
		else
			return 1;
	}       

	inline void clear(void)	
		// ZERO-OUT ALL ELEMENTS
	{
		while ( this->is_not_empty() )
		{
			this->erase_last_one();
		}
	}       
	inline void reset(void)	
		// INITIALIZE THE LIST, DELETES 'PTR_NEXT'
	{
		this->clear();

		delete ptr_next;
		ptr_next = NULL;

		tot_dim = MAX;
		last_one = -1;
	}       
				
	inline T &operator[] (int i)
	{
		if (i<0)
		{
			cout << "Negative index in the array...!\n";
			exit(1);
		}

		infobjarray<T,MAX>* ptr_current = this;

		while ( i >= MAX )
		{
			if (ptr_current->ptr_next == NULL)
			{
				ptr_current->ptr_next = new infobjarray<T,MAX>;
				tot_dim += MAX;
			}
			ptr_current = ptr_current->ptr_next;
			i -= MAX;
		}

		return ptr_current->data[i];
	}       

	infobjarray<T,MAX>* ptr_next;
	T* data;

	int last_one;
	int tot_dim;
};

#endif // INFOBJARRAY_H
