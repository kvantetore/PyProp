#ifndef PYPROP_DYNAMIC_CAST
#define PYPROP_DYNAMIC_CAST

#include <boost/shared_ptr.hpp>

template<class DestType, class SourceType>
boost::shared_ptr<DestType> pyprop_dynamic_cast(boost::shared_ptr<SourceType> sourcePtr)
{
#ifdef BZ_DEBUG
	shared_ptr<DestType> destPtr = boost::dynamic_pointer_cast<DestType>(sourcePtr);
	if (destPtr == 0)
	{
		cout 
			<< "WARNING: Could not cast " << typeid(*sourcePtr).name() 
			<< " to " << typeid(DestType).name() << " safely" << endl;
		destPtr = boost::static_pointer_cast<DestType>(sourcePtr);
	}
	return destPtr;
#else
	return boost::static_pointer_cast<DestType>(sourcePtr);
#endif
}

#endif
