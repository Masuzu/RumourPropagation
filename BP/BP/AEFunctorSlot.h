#pragma once

#include <memory>

//==================================================//
// Helper classes to wrap functions / class methods //
//==================================================//

class AEFunctorSlotBase
{
public:
	virtual ~AEFunctorSlotBase() {}
	virtual void operator()(void *, void *, void *)	const	{}
	virtual void operator()(void *, void *)	const	{}
	virtual void operator()(void *)	const = 0;
	virtual void operator()()	const = 0;
};

template<typename T, typename U>
class AEMemberFunctionSlot1 : public AEFunctorSlotBase
{
public:
	//! Member function slot type.
	typedef void(__stdcall T::*MemberFunctionType)(U);

	AEMemberFunctionSlot1(MemberFunctionType func, T* obj, U arg) :
		d_function(func),
		d_object(obj),
		d_arg(arg)
	{}

	AEMemberFunctionSlot1(MemberFunctionType func, T* obj) :
		d_function(func),
		d_object(obj)
	{}

	virtual void operator()()	const
	{
		return (d_object->*d_function)(d_arg);
	}

	virtual void operator()(void *arg)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg));
	}

private:
	MemberFunctionType d_function;
	T* d_object;
	U d_arg;
};

template<typename T, typename U, typename V>
class AEMemberFunctionSlot2 : public AEFunctorSlotBase
{
public:
	//! Member function slot type.
	typedef void(__stdcall T::*MemberFunctionType)(U, V);

	AEMemberFunctionSlot2(MemberFunctionType func, T* obj, U arg1, V arg2) :
		d_function(func),
		d_object(obj),
		d_arg1(arg1),
		d_arg2(arg2)
	{}

	AEMemberFunctionSlot2(MemberFunctionType func, T* obj, U arg1) :
		d_function(func),
		d_object(obj),
		d_arg1(arg1)
	{}

	AEMemberFunctionSlot2(MemberFunctionType func, T* obj) :
		d_function(func),
		d_object(obj)
	{}

	virtual void operator()()	const
	{
		return (d_object->*d_function)(d_arg1, d_arg2);
	}

	virtual void operator()(void *arg)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg), d_arg2);
	}

	virtual void operator()(void *arg1, void *arg2)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg1), reinterpret_cast<V>(arg2));
	}

private:
	MemberFunctionType d_function;
	T* d_object;
	U d_arg1;
	V d_arg2;
};

template<typename T, typename U, typename V, typename W>
class AEMemberFunctionSlot3 : public AEFunctorSlotBase
{
public:
	//! Member function slot type.
	typedef void(__stdcall T::*MemberFunctionType)(U, V, W);

	AEMemberFunctionSlot3(MemberFunctionType func, T* obj, U arg1, V arg2, W arg3) :
		d_function(func),
		d_object(obj),
		d_arg1(arg1),
		d_arg2(arg2),
		d_arg3(arg3)
	{}

	AEMemberFunctionSlot3(MemberFunctionType func, T* obj, U arg1, V arg2) :
		d_function(func),
		d_object(obj),
		d_arg1(arg1),
		d_arg2(arg2)
	{}

	AEMemberFunctionSlot3(MemberFunctionType func, T* obj, U arg1) :
		d_function(func),
		d_object(obj),
		d_arg1(arg1)
	{}

	AEMemberFunctionSlot3(MemberFunctionType func, T* obj) :
		d_function(func),
		d_object(obj)
	{}

	virtual void operator()()	const
	{
		return (d_object->*d_function)(d_arg1, d_arg2, d_arg3);
	}

	virtual void operator()(void *arg)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg), d_arg2, d_arg3);
	}

	virtual void operator()(void *arg1, void *arg2)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg1), reinterpret_cast<V>(arg2), d_arg3);
	}

	virtual void operator()(void *arg1, void *arg2, void *arg3)	const
	{
		return (d_object->*d_function)(reinterpret_cast<U>(arg1), reinterpret_cast<V>(arg2), reinterpret_cast<W>(arg3));
	}

private:
	MemberFunctionType d_function;
	T* d_object;
	U d_arg1;
	V d_arg2;
	W d_arg3;
};

template<typename T, typename U>
class AEFunctionPointerSlot1 : public AEFunctorSlotBase
{
public:
	AEFunctionPointerSlot1(T* function, U arg) :
		d_function(function),
		d_arg(arg)
	{}

	virtual void operator()()	const
	{
		return (*d_function)(d_arg);
	}

	virtual void operator()(void *arg)	const
	{
		return (*d_function)(reinterpret_cast<U>(arg));
	}

private:
	T* d_function;
	U d_arg;
};

class AEFunctor
{
public:
	AEFunctor()
	{
	}

	AEFunctor(const AEFunctor &foo)
	{
		d_functor_impl = foo.d_functor_impl;
	}

	void operator()()	const
	{
		return (*d_functor_impl)();
	}

	template<typename U>
	void operator()(U arg)	const
	{
		return (*d_functor_impl)(arg);
	}

	template<typename U, typename V>
	void operator()(U arg1, V arg2)	const
	{
		return (*d_functor_impl)(arg1, arg2);
	}

	template<typename U, typename V, typename W>
	void operator()(U arg1, V arg2, W arg3)	const
	{
		return (*d_functor_impl)(arg1, arg2, arg3);
	}

	template<typename T, typename U, typename V>
	AEFunctor(void (__stdcall T::*function)(U, V), T* obj, U arg1, V arg2) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot2<T, U, V>(function, obj, arg1, arg2)) )
	{}

	template<typename T, typename U, typename V>
	AEFunctor(void (__stdcall T::*function)(U, V), T* obj, U arg) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot2<T, U, V>(function, obj, arg)) )
	{}

	template<typename T, typename U, typename V>
	AEFunctor(void (__stdcall T::*function)(U, V), T* obj) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot2<T, U, V>(function, obj)) )
	{}

	template<typename T, typename U, typename V, typename W>
	AEFunctor(void (__stdcall T::*function)(U, V, W), T* obj, U arg1, V arg2, W arg3) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot3<T, U, V, W>(function, obj, arg1, arg2, arg3)) )
	{}

	template<typename T, typename U, typename V, typename W>
	AEFunctor(void (__stdcall T::*function)(U, V, W), T* obj, U arg1, V arg2) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot3<T, U, V, W>(function, obj, arg1, arg2)) )
	{}

	template<typename T, typename U, typename V, typename W>
	AEFunctor(void (__stdcall T::*function)(U, V, W), T* obj, U arg1) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot3<T, U, V, W>(function, obj, arg1)) )
	{}

	template<typename T, typename U, typename V, typename W>
	AEFunctor(void (__stdcall T::*function)(U, V, W), T* obj) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot3<T, U, V, W>(function, obj)) )
	{}

	template<typename T, typename U>
	AEFunctor(void (__stdcall T::*function)(U), T* obj, U arg) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot1<T, U>(function, obj, arg)) )
	{}

	template<typename T, typename U>
	AEFunctor(void (__stdcall T::*function)(U), T* obj) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEMemberFunctionSlot1<T, U>(function, obj)) )
	{}

	template<typename T, typename U>
	AEFunctor(T* function, U arg) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEFunctionPointerSlot1<T, U>(function, arg)) )
	{}

	template<typename T, typename U>
	AEFunctor(T* function) :
		d_functor_impl( std::shared_ptr<AEFunctorSlotBase>(new AEFunctionPointerSlot1<T, U>(function)) )
	{}

	void __stdcall Reset()
	{
		d_functor_impl.reset();
	}

private:
	//! Points to the internal functor object to which we are bound.
	std::shared_ptr<AEFunctorSlotBase> d_functor_impl;
};
