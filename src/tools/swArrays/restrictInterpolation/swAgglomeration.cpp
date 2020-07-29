#include "swAgglomeration.hpp"

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

namespace UNAP
{
	//- restrict
	bool* swRestInterMap::restFirstUse_ = NULL;
	restStruct* swRestInterMap::restStructLevels_ = NULL;

	//- face restrict
	bool* swRestInterMap::faceRestFirstUse_ = NULL;
	restStruct* swRestInterMap::faceRestStructLevels_ = NULL;

	//- interpolate
	bool* swRestInterMap::interFirstUse_ = NULL;
	interStruct* swRestInterMap::interStructLevels_ = NULL;

	//- agglomerate matrix upper
	bool* swRestInterMap::aggMatrixUpperFirstUse_ = NULL;
	aggMatrixUpperStruct* swRestInterMap::aggMatrixUpperStructLevels_ = NULL;


	label swRestInterMap::bandSize_ = 6250;
	label swRestInterMap::minCellsUsingSW_ = 10000;
}


void UNAP::swRestInterMap::initRestInterSize()
{
	const label levelSize = aggl_.size();

	if(!restFirstUse_)
	{
		restFirstUse_ = new bool [levelSize];
		restStructLevels_ = new restStruct [levelSize];

		forAll(i, levelSize)
		{
			restFirstUse_[i] = true;
		}
	}

	if(!faceRestFirstUse_)
	{
		faceRestFirstUse_ = new bool [levelSize];
		faceRestStructLevels_ = new restStruct [levelSize];

		forAll(i, levelSize)
		{
			faceRestFirstUse_[i] = true;
		}
	}


	if(!interFirstUse_)
	{
		interFirstUse_ = new bool [levelSize];
		interStructLevels_ = new interStruct [levelSize];

		forAll(i, levelSize)
		{
			interFirstUse_[i]  = true;
		}
	}

	if(!aggMatrixUpperFirstUse_)
	{
		aggMatrixUpperFirstUse_ = new bool [levelSize];
		aggMatrixUpperStructLevels_ = new aggMatrixUpperStruct [levelSize];

		forAll(i, levelSize)
		{
			aggMatrixUpperFirstUse_[i]  = true;
		}
	}
}


void UNAP::swRestInterMap::initRestStruct
(
	Vector<scalar>& cf,
	const Vector<scalar>& ff,
	const label fineLevelIndex
)
{
#ifdef SWTIMER
	swTimer::startTimer("restrictField init");
#endif
	if(restFirstUse_[fineLevelIndex])
	{
		const swInt* restrictMap  = aggl_.restrictAddressing(fineLevelIndex).begin();
		const swInt fSize = ff.size();
		const swInt cSize = cf.size();
		swInt slaveCycles = cSize / (bandSize_*64);
	    if(cSize % (bandSize_*64)) ++slaveCycles;
	    swInt slaveCores  = 64 * slaveCycles;

	    swInt remainder = cSize % slaveCores;
	    swInt lenShort  = cSize / slaveCores;
	    swInt lenLong   = lenShort + 1;

	    swInt fRemainder = fSize % slaveCores;
	    swInt fLenShort  = fSize / slaveCores;
	    swInt fLenLong   = fLenShort + 1;

	    //- allocate range
	    swInt** range = new swInt*[slaveCores];
	    forAll(i, slaveCores)
	    {
	        range[i] = new swInt[4];
	    }

	    forAll(i, slaveCores)
	    {
	        if (i < remainder)
	        {
	            range[i][0] = i * lenLong;
	            range[i][1] = range[i][0] + lenLong - 1;
	        }
	        else
	        {
	            range[i][0] = i * lenShort + remainder;
	            range[i][1] = range[i][0] + lenShort - 1;
	        }

	        //- range of fine data estimated
	        if (i < fRemainder)
	        {
	            range[i][2] = i * fLenLong;
	            range[i][3] = range[i][2] + fLenLong - 1;
	        }
	        else
	        {
	            range[i][2] = i * fLenShort + fRemainder;
	            range[i][3] = range[i][2] + fLenShort - 1;
	        }
	    }

	    //- check how many points will not be computed in slave cores
	    forAll(i, fSize)
	    {
	    	swInt cPos = restrictMap[i];
	        swInt cSlavePos = -1;
	        if(cPos < remainder*lenLong)
	        {
	            cSlavePos = cPos / lenLong;
	        }
	        else
	        {
	            cSlavePos = (cPos - remainder*lenLong) / lenShort + remainder;
	        }

	        if(i < range[cSlavePos][2])
	        {
	            range[cSlavePos][2] = range[cSlavePos][2] < i? range[cSlavePos][2] : i;
	        }
	        else if(i > range[cSlavePos][3])
	        {
	            range[cSlavePos][3] = range[cSlavePos][3] > i? range[cSlavePos][3] : i;
	        }
	    }

	    restStructLevels_[fineLevelIndex].mapPtr = restrictMap;
	    restStructLevels_[fineLevelIndex].localStartEnd = range;
	    restStructLevels_[fineLevelIndex].slaveCycles   = slaveCycles;

		restFirstUse_[fineLevelIndex] = false;
	}

	restStructLevels_[fineLevelIndex].fPtr = ff.begin();
	restStructLevels_[fineLevelIndex].cPtr = cf.begin();
#ifdef SWTIMER
	swTimer::endTimer("restrictField init");
#endif
}


void UNAP::swRestInterMap::initInterStruct
(
	Vector<scalar>& ff,
    const Vector<scalar>& cf,
    const label levelIndex
)
{
#ifdef SWTIMER
	swTimer::startTimer("prolongField init");
#endif
	if(interFirstUse_[levelIndex])
	{
		const swInt* restrictMap  = aggl_.restrictAddressing(levelIndex).begin();
		const swInt fSize = ff.size();
		const swInt cSize = cf.size();

		//- translate restrictMap
		swInt* interpolateMap = new swInt[fSize];
		swInt* interMapOffset = new swInt[cSize+1];

		swInt* offsetTemp = new swInt[cSize];

		forAll(i, cSize)
		{
			offsetTemp[i] = 0;
		}

		forAll(i, fSize)
		{
			swInt cPos = restrictMap[i];
			offsetTemp[cPos]++;
		}

		interMapOffset[0] = 0;
		forAll(i, cSize)
		{
			interMapOffset[i+1] = interMapOffset[i] + offsetTemp[i];
			offsetTemp[i] = 0;
		}

		forAll(i, fSize)
		{
			swInt cPos = restrictMap[i];
			interpolateMap[interMapOffset[cPos]+offsetTemp[cPos]] = i;
			offsetTemp[cPos]++;
		}

		delete []offsetTemp;

		swInt slaveCycles = fSize / (bandSize_*64);
		if(fSize % (bandSize_*64)) ++slaveCycles;
		swInt slaveCores  = 64 * slaveCycles;

	    swInt remainder = fSize % slaveCores;
	    swInt lenShort  = fSize / slaveCores;
	    swInt lenLong   = lenShort + 1;

	    swInt cRemainder = cSize % slaveCores;
	    swInt cLenShort  = cSize / slaveCores;
	    swInt cLenLong   = cLenShort + 1;

	    //- allocate range
	    swInt** range = new swInt*[slaveCores];
	    forAll(i, slaveCores)
	    {
	        range[i] = new swInt[4];
	    }

	    forAll(i, slaveCores)
	    {
	        //- range of fine data
	        if (i < remainder)
	        {
	            range[i][0] = i * lenLong;
	            range[i][1] = range[i][0] + lenLong - 1;
	        }
	        else
	        {
	            range[i][0] = i * lenShort + remainder;
	            range[i][1] = range[i][0] + lenShort - 1;
	        }

	        //- range of coarse data estimated
	        if(i < cRemainder)
	        {
	        	range[i][2] = i * cLenLong;
	        	range[i][3] = range[i][2] + cLenLong - 1;
	        }
	        else
	        {
	        	range[i][2] = i * cLenShort + cRemainder;
	            range[i][3] = range[i][2] + cLenShort - 1;
	        }
	    }



	    //- check how many points will not be computed in slave cores
	    forAll(i, cSize)
	    {
	    	swInt fCells = interMapOffset[i+1] - interMapOffset[i];
	    	swInt fMax = -1;
	    	swInt fMin = 1e+8;

	    	forAll(j, fCells)
	    	{
	    		swInt fPos = interpolateMap[interMapOffset[i]+j];
	    		fMax = fMax < fPos? fPos : fMax;
	    		fMin = fMin > fPos? fPos : fMin;
	    	}

	    	swInt fMaxSlavePos = -1, fMinSlavePos = -1;
	    	if(fMax < remainder*lenLong)
	    	{
	    		fMaxSlavePos = fMax / lenLong;
	    	}
	    	else
	    	{
	    		fMaxSlavePos = (fMax - remainder*lenLong) / lenShort + remainder;
	    	}

	    	if(fMin < remainder*lenLong)
	    	{
	    		fMinSlavePos = fMin / lenLong;
	    	}
	    	else
	    	{
	    		fMinSlavePos = (fMin - remainder*lenLong) / lenShort + remainder;
	    	}


	    	if(i < range[fMaxSlavePos][2])
	    	{
	    		range[fMaxSlavePos][2] = range[fMaxSlavePos][2] < i? range[fMaxSlavePos][2] : i;
	    	}
	    	else if(i > range[fMaxSlavePos][3])
	    	{
	    		range[fMaxSlavePos][3] = range[fMaxSlavePos][3] > i? range[fMaxSlavePos][3] : i;
	    	}


	    	if(i < range[fMinSlavePos][2])
	    	{
	    		range[fMinSlavePos][2] = range[fMinSlavePos][2] < i? range[fMinSlavePos][2] : i;
	    	}
	    	else if(i > range[fMinSlavePos][3])
	    	{
	    		range[fMinSlavePos][3] = range[fMinSlavePos][3] > i? range[fMinSlavePos][3] : i;
	    	}
	    }

	    interStructLevels_[levelIndex].mapPtr = interpolateMap;
    	interStructLevels_[levelIndex].offsetMapPtr = interMapOffset;
    	interStructLevels_[levelIndex].localStartEnd = range;
    	interStructLevels_[levelIndex].slaveCycles   = slaveCycles;

    	interFirstUse_[levelIndex] = false;
	}

	interStructLevels_[levelIndex].fPtr = ff.begin();
	interStructLevels_[levelIndex].cPtr = cf.begin();

#ifdef SWTIMER
	swTimer::endTimer("prolongField init");
#endif
}

namespace UNAP
{
template<>
void matrix::agglomeration::restrictField
(
    Vector<scalar>& cf,
    const Vector<scalar>& ff,
    const label fineLevelIndex
) const
{
    label len = ff.size();
    cf.SET_zero();

#ifdef SWTIMER
	swTimer::startTimer("restrictField");
#endif
    if(len > swRestInterMap::minCellsUsingSW_)
    {
    	swRestInterMap restMap(*this);
    	restMap.initRestInterSize();
    	restMap.initRestStruct(cf, ff, fineLevelIndex);
    	restrictData_host(&swRestInterMap::restStructLevels_[fineLevelIndex]);
    }
    else
    {
    	const labelVector& fineToCoarse = restrictAddressing(fineLevelIndex);
    	forAll(i, len)
	    {
	        cf[fineToCoarse[i]] += ff[i];
	    }
    }
#ifdef SWTIMER
	swTimer::endTimer("restrictField");
#endif
}


template<>
void matrix::agglomeration::prolongField
(
    Vector<scalar>& ff,
    const Vector<scalar>& cf,
    const label coarseLevelIndex
) const
{
    label len = ff.size();
    const labelVector &fineToCoarse = restrictAddressing(coarseLevelIndex);

#ifdef SWTIMER
	swTimer::startTimer("prolongField");
#endif
    if(len > swRestInterMap::minCellsUsingSW_)
    {
    	swRestInterMap interMap(*this);
    	interMap.initRestInterSize();
    	interMap.initInterStruct(ff, cf, coarseLevelIndex);
    	interpolateData_host(&swRestInterMap::interStructLevels_[coarseLevelIndex]);
    }
    else
    {
    	forAll(i, len)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
#ifdef SWTIMER
	swTimer::endTimer("prolongField");
#endif
}
}
