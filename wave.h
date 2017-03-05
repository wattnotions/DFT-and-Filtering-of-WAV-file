// WAVE file header format
struct HEADER {
    unsigned char riff[4];                      // RIFF string
    unsigned int overallSize   ;               // overall size of file in bytes
    unsigned char wave[4];                      // WAVE string
    unsigned char fmtChunkMarker[4];          // fmt string with trailing null char
    unsigned int lengthOfFmt;                 // length of the format data
    unsigned int formatType;                   // format type. 1-PCM, 3- IEEE float, 6 - 8bit A law, 7 - 8bit mu law
    unsigned int channels;                      // no.of channels
    unsigned int sampleRate;                   // sampling rate (blocks per second)
    unsigned int byterate;                      // SampleRate * NumChannels * BitsPerSample/8
    unsigned int blockAlign;                   // NumChannels * BitsPerSample/8
    unsigned int bitsPerSample;               // bits per sample, 8- 8bits, 16- 16 bits etc
    unsigned char dataChunkHeader [4];        // DATA string or FLLR string
    unsigned int dataSize;                     // NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
};

