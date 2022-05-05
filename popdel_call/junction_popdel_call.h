#ifndef JUNCTION_POPDEL_CALL_H_
#define JUNCTION_POPDEL_CALL_H_

using namespace seqan;
// =======================================================================================
// Enum Struct SVType
// =======================================================================================
// Enum for variant types
enum struct SVType: uint8_t
{
    UNKNOWN = 0,
    DEL = 1,
    DUP = 2,
    INV = 3,
    TRL = 4
};
std::ostream & operator << (std::ostream & out, const SVType & t)
{
    if (t == SVType::DEL )
        out << "DEL";
    else if (t == SVType::DUP )
        out << "DUP";
    else if (t == SVType::INV )
        out << "INV";
    else if (t == SVType::TRL )
        out << "BND";
    else
        out << "UNKNOWN";
    return out;
}
// =======================================================================================
// Function orientationToSVType
// =======================================================================================
// Return the SVType infered from the orientation.
inline SVType orientationToSVType(const Orientation & o)
{
    if (o == Orientation::FR)
        return SVType::DEL;
    if (o == Orientation::RF)
        return SVType::DUP;
    if (o == Orientation::FF || o == Orientation::RR)
        return SVType::INV;
}
// =======================================================================================
// Function svTypeToString
// =======================================================================================
// Return the SVType as a CharString
inline CharString svTypeToString(const SVType & t)
{
    if (t == SVType::DEL )
        return "DEL";
    else if (t == SVType::DUP )
        return "DUP";
    else if (t == SVType::INV )
        return "INV";
    else if (t == SVType::TRL )
        return "TRL";
    else
        return "UNKNOWN";
}
// =======================================================================================
// Struct JunctionCall
// =======================================================================================
// Hold information on a junction call
struct JunctionCall
{
    unsigned windowPosition;     // Starting position of the window used to identify the junction.
    String<Pair<unsigned> >perSampleFRPositions;    // i1: Local; i2: Mate chromsome. Used for translocations.
    String<Pair<unsigned> >perSampleRFPositions;    // i1: Local; i2: Mate chromsome. Used for translocations.
    String<Pair<unsigned> >perSampleFFPositions;    // i1: Local; i2: Mate chromsome. Used for translocations.
    String<Pair<unsigned> >perSampleRRPositions;    // i1: Local; i2: Mate chromsome. Used for translocations.
    unsigned position;           // Position on the chromsome
    int mateChromosome;           // Only for translocations: Chromsome of the mate
    unsigned matePosition;       // Estimated position of the associated junction. 0 if unwknown
    Orientation orientation;     // Orientation of the supporting read pairs.
    SVType svtype;               // Predicted type of the variant belonging to the junction
    String<Tuple<double, 5> > perSampleTranslocSupport;// Number of support read pairs for each sample for translocations.
                                             // i[0]:FR, 1:rf, 2:ff, 3:rr, 4:Non-supporting. If merged: sum of windows
    String<double> perSampleSupport;  // Number of supporting read pairs for each sample. If merged: Avg. per window.
    double totalSupport;       // Number off all supporting reads pairs across all samples
    String<Triple<double> > gtLikelihoods; // GT likelihhods of each sample. Order: HomRef, Het, HomDel.
    double frequency;             // Allele-frequency of the deletion across all samples.
    String<float> perSampleCoverageChange; // Relative change in coverage compared to previous window per sample. 1.0 indicates no change.
    //Insert sizes ?
     //Clipping behavior ?

    private:
        void init(unsigned wpos,
                  String<Pair<unsigned> > & psfrp,
                  String<Pair<unsigned> > & psrfp,
                  String<Pair<unsigned> > & psffp,
                  String<Pair<unsigned> > & psrrp,
                  unsigned pos,
                  int mChrom,
                  unsigned mPos,
                  double f,
                  Orientation o,
                  SVType t,
                  String<Tuple<double, 5> > & psts,
                  String<double> & pss,
                  double ts,
                  String<Triple<double> > & gtl,
                  String<float> & pscc)
        {
            windowPosition = wpos;
            perSampleFRPositions = psfrp;
            perSampleRFPositions = psrfp;
            perSampleFFPositions = psffp;
            perSampleRRPositions = psrrp;
            position = pos;
            mateChromosome = mChrom;
            matePosition = mPos;
            frequency = f;
            orientation = o;
            svtype = t;
            perSampleTranslocSupport = psts;
            perSampleSupport = pss;
            totalSupport = ts;
            gtLikelihoods = gtl;
            perSampleCoverageChange = pscc;
        }
    public:
        JunctionCall()
        {
            String<Pair<unsigned> > psfrp;
            String<Pair<unsigned> > psrfp;
            String<Pair<unsigned> > psffp;
            String<Pair<unsigned> > psrrp;
            String<Tuple<double, 5> > psts;
            String<double> pss;
            String<Triple<double> > gtl;
            String<float> pscc;
            resize(psfrp, 0);
            resize(psfrp, 0);
            resize(psts, 0);
            resize(pss, 0);
            resize(gtl, 0);
            resize(pscc, 0);
            init(0, psfrp, psrfp, psffp, psrrp, 0, -1, 0, 0.0, Orientation::FR, SVType::UNKNOWN, psts, pss, 0, gtl, pscc);
        }
        JunctionCall(unsigned wpos)
        {
            String<Pair<unsigned> > psfrp;
            String<Pair<unsigned> > psrfp;
            String<Pair<unsigned> > psffp;
            String<Pair<unsigned> > psrrp;
            String<Tuple<double, 5> > psts;
            String<double> pss;
            String<Triple<double> > gtl;
            String<float> pscc;
            resize(psfrp, 0);
            resize(psrfp, 0);
            resize(psts, 0);
            resize(pss, 0);
            resize(gtl, 0);
            resize(pscc, 0);
            init(wpos, psfrp, psrfp, psffp, psrrp, 0, -1, 0, 0.0, Orientation::FR, SVType::UNKNOWN, psts, pss, 0, gtl, pscc);
        }
        JunctionCall(unsigned wpos,
                     String<Pair<unsigned> > & psfrp,
                     String<Pair<unsigned> > & psrfp,
                     String<Pair<unsigned> > & psffp,
                     String<Pair<unsigned> > & psrrp,
                     unsigned pos,
                     unsigned mPos,
                     double f,
                     Orientation o,
                     SVType t,
                     String<Tuple<double, 5> > & psts,
                     String<double> & pss,
                     double ts,
                     String<Triple<double> > & gtl,
                     String<float> & pscc)
        {
            init(wpos, psfrp, psrfp, psffp, psrrp, pos, -1, mPos, f, o, t, psts, pss, ts, gtl, pscc);
        }
        JunctionCall(unsigned wpos,
                     String<Pair<unsigned> > & psfrp,
                     String<Pair<unsigned> > & psrfp,
                     String<Pair<unsigned> > & psffp,
                     String<Pair<unsigned> > & psrrp,
                     unsigned pos,
                     unsigned mPos,
                     double f,
                     Orientation o,
                     String<Tuple<double, 5> > & psts,
                     String<double> & pss,
                     double ts,
                     String<Triple<double> > & gtl,
                     String<float> & pscc)
        {
            init(wpos, psfrp, psrfp, psffp, psrrp, pos, -1, mPos, f, o, orientationToSVType(o), psts, pss, ts, gtl, pscc);
        }
        JunctionCall(unsigned wpos,
                    String<Pair<unsigned> > & psfrp,
                    String<Pair<unsigned> > & psrfp,
                    String<Pair<unsigned> > & psffp,
                    String<Pair<unsigned> > & psrrp,
                    unsigned pos,
                    int mChrom,
                    unsigned mPos,
                    double f,
                    Orientation o,
                    String<Tuple<double, 5> > & psts,
                    String<double> & pss,
                    double ts,
                    String<Triple<double> > & gtl,
                    String<float> & pscc)
        {
            init(wpos, psfrp, psrfp, psffp, psrrp, pos, mChrom, mPos, f, o, orientationToSVType(o), psts, pss, ts, gtl, pscc);
        }
        // =======================================================================================
        // Function reset()
        // =======================================================================================
        // Resets the JunctionCall
        inline void reset(void)         // Maybe also resize all per sample info to 0 to save memory?
        {
            windowPosition = 0;
            clear(perSampleFRPositions);
            clear(perSampleRFPositions);
            clear(perSampleFFPositions);
            clear(perSampleRRPositions);
            position = 0;
            mateChromosome = -1;
            matePosition = 0;
            frequency = 0.0;
            orientation = Orientation::FR;
            svtype = SVType::UNKNOWN;
            clear(perSampleTranslocSupport);
            clear(perSampleSupport);
            totalSupport = 0;
            clear(gtLikelihoods);
            clear(perSampleCoverageChange);
        }
        // =======================================================================================
        // Function print()
        // =======================================================================================
        // Prints some information about the call. Only for debugging and testing purposes.
        inline void print(void) const
        {
            std::cout << "WPos: " << windowPosition <<"\tPosition" << ": " << position;
            if (svtype != SVType::TRL)
            {
                std::cout << "-" << matePosition << " (" << matePosition - position << "bp); Orientation/SVType: " 
                          << orientation << "/" << svtype;

                unsigned c = 0;
                for (auto it = begin(perSampleSupport); it != end(perSampleSupport); ++it)
                    if (*it > 0)
                        ++c;

                std::cout << "; Avg. Support: " << totalSupport << " pairs from "  << c << " sample(s):\t";
                for (auto it = begin(perSampleSupport); it != end(perSampleSupport); ++it)
                    std::cout << *it << "\t";

                std::cout << "\tGT likelihoods\t";
                for (auto it = begin(gtLikelihoods); it != end(gtLikelihoods); ++it)
                {
                    std::cout << it->i1 << ":" << it->i2 << ":" << it->i3;
                    if (it != end(gtLikelihoods) - 1)
                        std::cout << "\t";
                }
                std::cout << std::endl;

            }
            else    //Translocation
            {
                std::cout << "+(" << mateChromosome << "):" << matePosition << "; Orientation/SVType: " 
                          << orientation << "/" << svtype;
            }
        }
};
// -----------------------------------------------------------------------------
// Function isValid
// -----------------------------------------------------------------------------
// Return true if a Junction Call is valid.
inline bool isValid(const JunctionCall & j)
{
    return j.totalSupport > 0u;
}
// -----------------------------------------------------------------------------
// Function invalidateJunctionCall
// -----------------------------------------------------------------------------
// Marks the JunctionCall as invalid.
inline void invalidateJunctionCall(JunctionCall & j)
{
    j.totalSupport = 0u;
}
// -----------------------------------------------------------------------------
// Function printJunctionCalls
// -----------------------------------------------------------------------------
// Prints the window wise junction calls.
// Mainly for testing purposes.
inline void printJunctionCalls(const String<JunctionCall> & calls)
{
    for (auto it = begin(calls); it != end(calls); ++it)
    {
        if (isValid(*it))
            it->print();
    }
}
inline void printJunctionCalls(const Tuple<String<JunctionCall>, 4> & t)
{
    for (unsigned i = 0; i < 4; ++i)
        printJunctionCalls(t[i]);
}
// -----------------------------------------------------------------------------
// Function sameOrientation
// -----------------------------------------------------------------------------
// Return true if the two JunctionCalls have the same orientation, false otherwise.
inline bool sameOrientation(const JunctionCall & l, const JunctionCall & r)
{
    return l.orientation == r.orientation;
}
// -----------------------------------------------------------------------------
// Function similarStart
// -----------------------------------------------------------------------------
// Return true if the two JunctionCalls have a similar position false otherwise.
inline bool similarPosition(const JunctionCall & l,
                            const JunctionCall & r,
                            const double & stddev,
                            const double & f = 3.0)
{
    return (std::max(l.position, r.position) - std::min(l.position, r.position)) < f * stddev;
}
// -----------------------------------------------------------------------------
// Function similarEnd
// -----------------------------------------------------------------------------
// Return true if the two JunctionCalls have a similar mate position false otherwise.
inline bool similarMatePosition(const JunctionCall & l,
                                const JunctionCall & r,
                                const double & stddev,
                                const double & f = 3.0)
{
    return std::max(l.matePosition, r.matePosition) - std::min(l.matePosition, r.matePosition) < f * stddev;
}
// -----------------------------------------------------------------------------
// Function similar
// -----------------------------------------------------------------------------
// Return true if the two JunctionCalls are similar enough to be considered the same.
inline bool similar(const JunctionCall & l, const JunctionCall & r, const double & stddev, const double & f = 4.0)
{
    if (!sameOrientation(l, r))
        return false;
    if (!similarPosition(l, r, stddev, f))
        return false;
    if (!similarMatePosition(l, r, stddev, f))
        return false;

    return true;
}
struct TranslocJunctions{};
inline bool similar(const JunctionCall & l,
                    const JunctionCall & r,
                    const double & stddev,
                    TranslocJunctions flag,
                    const double & f = 4.0)
{
    (void) flag;
    if (l.mateChromosome != r.mateChromosome)
        return false;
    if (!similarPosition(l, r, stddev, f))
        return false;
    if (!similarMatePosition(l, r, stddev, f))
        return false;

    return true;
}
// -----------------------------------------------------------------------------
// operator< overload
// -----------------------------------------------------------------------------
// Define the < operator for JunctionCall objects.
inline bool operator<(const JunctionCall & l, const JunctionCall & r)
{
    if (l.position < r.position)
        return true;
    else if (l.position > r.position)
        return false;
    else if (l.matePosition < r.matePosition)
        return true;
    else if (l.matePosition > r.matePosition)
        return false;
    else if (l.windowPosition < r.windowPosition) // Might not be necessary later, but preserves the order of the win
        return true;
    else
        return false;
}
// -----------------------------------------------------------------------------
// translocOrder
// -----------------------------------------------------------------------------
// Return true if the ttrranslocation call defined in l comes before r when ordering window wise calls
inline bool translocOrder(const JunctionCall & l, const JunctionCall & r)
{
    if (l.mateChromosome != r.mateChromosome)
        return l.mateChromosome < r.mateChromosome;
    else if (l.position != r.position)
        return l.position < r.position;
    else if (l.matePosition != r.matePosition)
        return l.matePosition < r.matePosition;
    else
        return l.windowPosition < r.windowPosition; // Might not be necessary later, but preserves the order of the win
}
// =======================================================================================
// Struct OrientationCounts
// =======================================================================================
// Hold counts of read pairs per orientation.
struct OrientationCounts
{
    unsigned frCount;
    unsigned rfCount;
    unsigned ffCount;
    unsigned rrCount;

    OrientationCounts()
    {
        frCount = 0;
        rfCount = 0;
        ffCount = 0;
        rrCount = 0;
    }
    OrientationCounts(unsigned fr, unsigned rf, unsigned ff, unsigned rr)
    {
        frCount = fr;
        rfCount = rf;
        ffCount = ff;
        rrCount = rr;
    }
    inline void reset(void)
    {
        frCount = rfCount = ffCount = rrCount = 0u;
    }
    // =======================================================================================
    // Function add()
    // =======================================================================================
    // Increment the counts for the given orientation by one.
    inline void add(const Orientation & o)
    {
        if (o == Orientation::FR)
            ++frCount;
        else if (o == Orientation::RF)
            ++rfCount;
        else if (o == Orientation::FF)
            ++ffCount;
        else // if (o == Orientation::RR)
            ++rrCount;
    }
    inline OrientationCounts& operator+=(const OrientationCounts & r)
    {
        this->frCount += r.frCount;
        this->rfCount += r.rfCount;
        this->ffCount += r.ffCount;
        this->rrCount += r.rrCount;
        return *this;
    }
    inline OrientationCounts& operator-=(const OrientationCounts & r)
    {
        SEQAN_ASSERT_GEQ(this->frCount, r.frCount);
        SEQAN_ASSERT_GEQ(this->rfCount, r.rfCount);
        SEQAN_ASSERT_GEQ(this->ffCount, r.ffCount);
        SEQAN_ASSERT_GEQ(this->rrCount, r.rrCount);
        this->frCount -= r.frCount;
        this->rfCount -= r.rfCount;
        this->ffCount -= r.ffCount;
        this->rrCount -= r.rrCount;
        return *this;
    }
};
inline OrientationCounts operator+(const OrientationCounts & l, const OrientationCounts & r)
{
    return OrientationCounts(l.frCount + r.frCount, l.rfCount + r.rfCount, l.ffCount + r.ffCount, l.rrCount +r.rrCount);
}
inline OrientationCounts operator-(const OrientationCounts & l, const OrientationCounts & r)
{
    SEQAN_ASSERT_GEQ(l.frCount, r.frCount);
    SEQAN_ASSERT_GEQ(l.rfCount, r.rfCount);
    SEQAN_ASSERT_GEQ(l.ffCount, r.ffCount);
    SEQAN_ASSERT_GEQ(l.rrCount, r.rrCount);
    return OrientationCounts(l.frCount - r.frCount, l.rfCount - r.rfCount, l.ffCount - r.ffCount, l.rrCount -r.rrCount);
}
inline bool operator==(const OrientationCounts & l, const OrientationCounts & r)
{
    return l.frCount == r.frCount && l.rfCount == r.rfCount && l.ffCount == r.ffCount && l.rrCount == r.rrCount;
}
inline unsigned sum(const OrientationCounts & o)
{
    return o.frCount + o.rfCount + o.ffCount + o.rrCount;
}
#endif /* JUNCTION_POPDEL_CALL_H_ */
