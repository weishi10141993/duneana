#ifndef PDSTriggerPrimitiveFinderTool_h
#define PDSTriggerPrimitiveFinderTool_h

#include <vector>
#include <iostream>

class PDSTriggerPrimitiveFinderTool {
 
public:
    struct Hit
    {
    Hit(int _channel, int _startTime, int _peakCharge, int _SADC, int _timeOverThreshold)
	: channel(_channel),
	startTime(_startTime),
	peakCharge(_peakCharge),
	SADC(_SADC),
	timeOverThreshold(_timeOverThreshold)
      {}
      int channel;
      int startTime;
      int peakCharge;
      int SADC;
      int timeOverThreshold;
    };

    virtual ~PDSTriggerPrimitiveFinderTool() =default;

    virtual std::vector<PDSTriggerPrimitiveFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples) = 0;
 
};

#endif // include guard
