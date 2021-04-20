#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int genChannelMap()
{
  string outputfile = "../../src/FstNoiseStudy/channelMap.txt";
  ofstream file_channelMap;
  file_channelMap.open(outputfile.c_str());

  if (!file_channelMap.is_open())
  {
    std::cout << "failed to open " << outputfile.c_str() << endl;
    return -1;
  }
  else
  {
    int readout_channelMap[128];
    int physics_channelMap[128];
    // readout channel vs. physics channel
    for(int i_ro = 0; i_ro < 128; ++i_ro)
    {
      int channel = 32*(i_ro%4) + 8*int(i_ro/4) - 31*int(i_ro/16);
      cout << "readout_channel = " << i_ro << ", physics_channel = " << channel << endl;
      // file_channelMap << i_ro << "    " << channel << endl;

      physics_channelMap[channel] = channel;
      readout_channelMap[channel] = i_ro;
    }

    // physics channel vs. readout channel
    for(int i_ch = 0; i_ch < 128; ++i_ch)
    {
      cout << "i_ch = " << i_ch << ", physics_channelMap = " << physics_channelMap[i_ch] << ", readout_channelMap = " << readout_channelMap[i_ch] << endl;
      file_channelMap << physics_channelMap[i_ch] << "    " << readout_channelMap[i_ch] << endl;
    }
  }

  file_channelMap.close();

  return 1;
}
