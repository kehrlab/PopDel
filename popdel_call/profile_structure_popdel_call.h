#ifndef PROFILE_STRUCTURE_POPDEL_H_
#define PROFILE_STRUCTURE_POPDEL_H_

#include <unordered_set>

#include <seqan/basic.h>
#include "../utils_popdel.h"
#include "../insert_histogram_popdel.h"

using namespace seqan;
// =======================================================================================
// Class ChromosomeProfileStartEntry
// =======================================================================================
// Stores the firstWindow and the insert size deviation of a specific insert.
struct ChromosomeProfileStartEntry
{
    __uint32 firstWindow;
    __int32  insertSizeDeviation;

    ChromosomeProfileStartEntry() :
        firstWindow(0), insertSizeDeviation(0){}

    ChromosomeProfileStartEntry(__uint32 start, __int32 deviation) :
        firstWindow(start), insertSizeDeviation(deviation){}
};
// =======================================================================================
// Operator== overload
// =======================================================================================
bool operator==(const ChromosomeProfileStartEntry & l, const ChromosomeProfileStartEntry & r)
{
    return (l.firstWindow == r.firstWindow && l.insertSizeDeviation == r.insertSizeDeviation);
}
bool operator==(const String<ChromosomeProfileStartEntry> & l, const String<ChromosomeProfileStartEntry> & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
bool operator==(const String<String<ChromosomeProfileStartEntry> > & l,
                const String<String<ChromosomeProfileStartEntry> > & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
// =======================================================================================
// Class StartEntrySubset()
// =======================================================================================
// Class representing one subset for the use in CyclicStartEntryTable.
struct StartEntrySubset
{
    String<ChromosomeProfileStartEntry>     subset;
    unsigned                                right;         // exclusiv right border.
    unsigned                                idOffset;
    bool                                    readable;      // True to indicate that the subset is open for reading.

    StartEntrySubset():
    right(0), idOffset(0), readable(false){}

    StartEntrySubset(unsigned idO, unsigned r) :
    right(r), idOffset(idO), readable(false){}
    // =======================================================================================
    // Function clear()
    // =======================================================================================
    //Clear the subset.
    inline void clear()
    {
        seqan::clear(subset);
        right = 0;
        idOffset = 0;
        readable = false;
    }
    // =======================================================================================
    // Function reset()
    // =======================================================================================
    // Clear the subset and prepare it to be used with the given offset and border.
    // Also make sure it is not used for reading.
    inline void reset(const unsigned & idO, const unsigned & rightBorder)
    {
        seqan::clear(subset);
        right = rightBorder;
        idOffset = idO;
        readable = false;
    }
    // =======================================================================================
    // Function permitReading()
    // =======================================================================================
    // Permit reading of the subset.
    inline void permitReading()
    {
        readable = true;
    }
    // =======================================================================================
    // Function forbidReading()
    // =======================================================================================
    // Forbid reading of the subset.
    inline void forbidReading()
    {
        readable = false;
    }
};
// =======================================================================================
// Operator== overloads for StartEntrySubset comparisons.
// =======================================================================================
bool operator==(const StartEntrySubset & l, const StartEntrySubset & r)
{
    if (length(l.subset) != length(r.subset))
        return false;
    if (l.right != r.right)
        return false;
    if (l.idOffset != r.idOffset)
        return false;
    if (l.readable != r.readable)
        return false;
    Iterator<const String<ChromosomeProfileStartEntry> >::Type lIt = begin(l.subset, Standard());
    Iterator<const String<ChromosomeProfileStartEntry> >::Type lEnd = end(l.subset, Standard());
    Iterator<const String<ChromosomeProfileStartEntry> >::Type rIt = begin(r.subset, Standard());
    while (lIt != lEnd)
    {
        if (!(*lIt == *rIt))
            return false;
        ++lIt;
        ++rIt;
    }
        return true;
}
bool operator==(const String<StartEntrySubset> & l, const String<StartEntrySubset> & r)
{
    if (length(l) != length(r))
        return false;
    Iterator<const String<StartEntrySubset> >::Type lIt  = begin(l, Standard());
    Iterator<const String<StartEntrySubset> >::Type lEnd = end(l, Standard());
    Iterator<const String<StartEntrySubset> >::Type rIt  = begin(r, Standard());
    while (lIt != lEnd)
    {
        if (!(*lIt == *rIt))
            return false;
        ++lIt;
        ++rIt;
    }
    return true;
}
// =======================================================================================
// Class CyclicStartEntryTable
// =======================================================================================
// Stores and manages the start entries of the chromosome profile.
// It is split into 3 sets, which are used in a round-robin-fashion.
// Entries are added to the subset with matching borders for the first window of the entrie.
// The values can be accessed using the []-operator with the ID (i.e. order of insertion) of the element.
// Only elements which have been added to one set-switch before may be accessed.
// Example: The CyclicStartEntryTable c(0, 10) allows entries with firstWindow between [0,10[ in the first set (sets[0]).
// After adding 5 elements which fulfill this condition to the set via c.add(firstWindow, x), one more element
// is added with c.add(13, x). As it lies outside of the first set's borders, it is added to sets[1] whose borders
// are set to [10,20[. Now (and not before the set-switch) the entries of firstSet may be read via c[id] (0<=id<5).
// After adding some 9 more elements to secondSet and one with c.add(23, x), which is assigned to sets[2],
// accessing sets[1] via c[id] (5<=id<15) becomes possible. If we now call c.add(35, x), which is added to sets[0],
// the first set only has one element: c[15] == (35, x). If we try c.add(55, x), the functions only shifs to the next
// set without adding the element, because it does not belong there. Calling c.add(55, x) again now add it to sets[2],
// with sets[1] being empty.
// Note: It is very imporant to keep the rythm of adding and reading from the set: 
// Keep adding while c.add() returns 0, then read using c.goNext() and c.nextRead unil the end of the subset is reached,
// then continue adding to the next set and so on.
struct CyclicStartEntryTable
{
    String<StartEntrySubset>                                    sets;
    Iterator<String<ChromosomeProfileStartEntry>, Rooted>::Type nextRead;
    unsigned                                                    readSet;        // Index of the set of nextRead.
    unsigned                                                    writeSet;
    unsigned                                                    insertionCount;  // Counts the number of add() calls.
    unsigned                                                    numWindows;      // Max number of windows per set.

    CyclicStartEntryTable(unsigned firstWin = 0, unsigned windows = 10000)
    {
        SEQAN_ASSERT_GT(windows, 0u);
        numWindows = windows;
        writeSet = 0;
        resize(sets, 3, Exact());
        sets[0].reset(0, firstWin + windows);
        sets[1].reset(0, sets[0].right + windows);
        sets[2].reset(0, sets[1].right + windows);
        for (unsigned i = 0; i < 3; ++i)
            reserve(sets[i].subset, numWindows);
        readSet = 2;
        nextRead = begin(sets[readSet].subset, Rooted());
        insertionCount = 0;
    }
    // =======================================================================================
    // Function getSetNumFromId()
    // =======================================================================================
    // Return the index in sets of the subset for the given ID.
    inline unsigned getSetNumFromId (const unsigned & id) const
    {
        if (writeSet == 0u)
        {
            if (sets[2].idOffset <= id)
                return 2;
            else
                return 1;
        }
        else if (writeSet == 1u)
        {
            if (sets[0].idOffset <= id)
                return 0;
            else
                return 2;
        }
        else //writeSet == 2
        {
            if (sets[1].idOffset <= id)
                return 1;
            else
                return 0;
        }
    }
    // =======================================================================================
    // Function correctId()
    // =======================================================================================
    // Correct the given ID, s.t. it points to the corresponding element in the readSet.
    // Return a reference to the set.
    inline StartEntrySubset& correctId(unsigned & id)
    {
        StartEntrySubset& set = sets[getSetNumFromId(id)];
        SEQAN_ASSERT_GEQ(id, set.idOffset);
        id -= set.idOffset;
        SEQAN_ASSERT(set.readable);
        return set;
    }
    // =======================================================================================
    // Function needsSwitch()
    // =======================================================================================
    // Return true if the given firstWindow requires a switch to a higher set, false otherwise.
    inline bool needsSwitch(const unsigned & firstWindow) const
    {
        if (firstWindow >= sets[writeSet].right)
            return true;
        else
            return false;
    }
    inline bool tooBigForNext(const unsigned & firstWindow) const
    {
        if (firstWindow >= sets[writeSet].right + numWindows)
            return true;
        else
            return false;
    }
    // =======================================================================================
    // Function switchWriteSet()
    // =======================================================================================
    // Switch between the three sets for writing. The borders and the bool indicating the active set. Also clear the
    // next set for writing.
    inline void switchWriteSet()
    {
        StartEntrySubset & previousSet = sets[writeSet];
        StartEntrySubset & nextSet = sets[(writeSet + 1) % 3];
        nextSet.reset(insertionCount, previousSet.right + numWindows);
        writeSet = (writeSet + 1) % 3;
        SEQAN_ASSERT_EQ(previousSet.readable, false);
        previousSet.permitReading();
    }
    // =======================================================================================
    // Function switchReadSet()
    // =======================================================================================
    // Set nextRead to the begin of the next set.
    // Return false if the nextSet is empty, true otherwise.
    inline bool switchReadSet()
    {
        readSet = (readSet + 1) % 3;
        SEQAN_ASSERT(sets[readSet].readable);
        nextRead = begin(sets[readSet].subset, Rooted());
        if (atEnd(nextRead))
            return false;
        else
            return true;
    }
    // =======================================================================================
    // Function backwardSwitchReadSet()
    // =======================================================================================
    // Set nextRead to the end-1 of the previous set.
    // Return false if the set is empty, true otherwise.
    inline bool backwardSwitchReadSet()
    {
        readSet = (readSet + 2) % 3;
        SEQAN_ASSERT(sets[readSet].readable);
        nextRead = end(sets[readSet].subset, Rooted());
        if (atBegin(nextRead))
        {
            return false;
        }
        else
        {
            --nextRead;
            return true;
        }
    }
    // =======================================================================================
    // Function switchBothSets()
    // =======================================================================================
    // Calls switchWriteSet() and switchReadSet().
    // This should alwasy result in nextRead being set to the beginning of the previous writeSet.
    // Eg. readSet == 1, writeSet == 2  -> readSet == 2; writeSet == 0;
    // Return false if the next set for reading is empty, true otherwise.
    inline bool switchBothSets(unsigned & currentRightBorder)
    {
        switchWriteSet();
        bool res = switchReadSet();
        currentRightBorder = sets[readSet].right;       // TODO: Move to another function to avoid overhead.
        return res;
    }
    // =======================================================================================
    // Function add()
    // =======================================================================================
    // Add the entry to the current writeSet and increase the insertionCount;
    // This function does not care if the currentSet is actually the correct writeSet!
    inline void add(const unsigned & firstWin, const int & deviation)
    {
        SEQAN_ASSERT_EQ(sets[writeSet].readable, false);    // Never insert into a readable subset.
        appendValue(sets[writeSet].subset, ChromosomeProfileStartEntry(firstWin, deviation));
        ++insertionCount;
    }
    // =======================================================================================
    // Operator[] overload
    // =======================================================================================
    // Return a reference to the i-th entry.
    // Reading only occurs from the currently INactive set. E.g. if firstSetActive == true, reading from secondSet.
    // There is no check, if the value at the desired position is valid or deprecated.
    inline ChromosomeProfileStartEntry & operator[](unsigned id)
    {
        StartEntrySubset& readSet = correctId(id);
        return(readSet.subset[id]);
    }
    // Move the iterator for reading forward by one.
    // Return false the iterator is at the end of the current subset.
    inline bool goNext()
    {
        SEQAN_ASSERT(!atEnd(nextRead));
        ++nextRead;
        if (atEnd(nextRead))
            return false;
        else
            return true;
    }
    // =======================================================================================
    // Function goPrevious()
    // =======================================================================================
    // Move the iterator for reading back by one and switch the set back by one if necessary.
    // Return true on succes and false if the set has been switch but the new set is empty.
    inline bool goPrevious()
    {
        SEQAN_ASSERT(sets[readSet].readable);
        if (atBegin(nextRead))
        {
            return backwardSwitchReadSet();
        }
        else
        {
            --nextRead;
            return true;
        }
    }
    // =======================================================================================
    // Function peekPreviousFirstWindow()
    // =======================================================================================
    // Return the firstWindow of the previous position without changing the iterator.
    // Return 0 if the previous position is invalid.
    inline unsigned peekPreviousFirstWindow() const
    {
        if (atBegin(nextRead))
        {
            unsigned prevSetNum = (readSet + 2) % 3;
            if (empty(sets[prevSetNum].subset) || sets[readSet].right < sets[prevSetNum].right)
                return 0;
            else
                return  back(sets[prevSetNum].subset).firstWindow;
        }
        else
        {
            return (nextRead - 1)->firstWindow;
        }
    }
    // =======================================================================================
    // Function getEntry()
    // =======================================================================================
    // Return a reference to the element pointed at by nextRead.
    inline const ChromosomeProfileStartEntry & getEntry() const
    {
        SEQAN_ASSERT(!atEnd(nextRead));
        return *nextRead;
    }
    // =======================================================================================
    // Function getEntryAt()
    // =======================================================================================
    // Return a reference to the element with the given uncorrected ID.
    inline const ChromosomeProfileStartEntry & getEntryAt(unsigned id) const
    {
        unsigned setNum = getSetNumFromId(id);
        SEQAN_ASSERT_LEQ(sets[setNum].idOffset, id);
        id -= sets[setNum].idOffset;
        return(sets[setNum].subset[id]);
    }
    // =======================================================================================
    // Function getDeviation()
    // =======================================================================================
    // Return the insert size deviation of nextRead.
    inline int getDeviation() const
    {
        return nextRead->insertSizeDeviation;
    }
    // =======================================================================================
    // Function getDeviationAt()
    // =======================================================================================
    // Return the insert size deviation of the element with the given ID.
    inline int getDeviationAt(unsigned id) const
    {
        return getEntryAt(id).insertSizeDeviation;
    }
    // =======================================================================================
    // Function getFirstWin()
    // =======================================================================================
    // Return firstWindow of nextRead.
    inline unsigned getFirstWin() const
    {
        return nextRead->firstWindow;
    }
    // =======================================================================================
    // Function getFirstWinAt()
    // =======================================================================================
    // Return firstWindow of the elment with the given ID.
    inline unsigned getFirstWinAt(unsigned id) const
    {
        return getEntryAt(id).firstWindow;
    }
    // =======================================================================================
    // Function getEntryId()
    // =======================================================================================
    // Return the uncorrected index (i.e. the ID) of the element pointed at by nextRead.
    inline unsigned getEntryId() const
    {
        return position(nextRead) + sets[readSet].idOffset;
    }
};
// =======================================================================================
// Operator== overloads
// =======================================================================================
bool operator==(const CyclicStartEntryTable & l, const CyclicStartEntryTable & r)
{
    if (!(l.sets == r.sets))
        return false;
    if (position(l.nextRead) != position(r.nextRead))
        return false;
    if (l.readSet != r.readSet)
        return false;
    if (l.writeSet != r.writeSet)
        return false;
    if (l.insertionCount != r.insertionCount)
        return false;
    if (l.numWindows != r.numWindows)
        return false;
    else
        return true;
}
bool operator==(const String<CyclicStartEntryTable> & l, const String<CyclicStartEntryTable> & r)
{
    if (length(l) != length(r))
        return false;
    bool res = true;
    for (unsigned i = 0; i < length(l); ++i)
    {
        res = l[i] == r[i];
        if (!res)
            return false;
    }
    return true;
}
// =======================================================================================
// Overload for Iterator functions atEnd() and atBegin()
// =======================================================================================
bool atEnd(const CyclicStartEntryTable & c)
{
    return atEnd(c.nextRead);
}
bool atBegin(const CyclicStartEntryTable & c)
{
    return atBegin(c.nextRead);
}
// =======================================================================================
// Function switchSets()
// =======================================================================================
// Wrapper for calling switchBothSets() on a String of CyclicStatEntryTable.
inline void switchSets(String<CyclicStartEntryTable> & s)
{
    unsigned tmp;
    for (Iterator<String<CyclicStartEntryTable> >::Type it = begin(s); it != end(s); ++it)
        it->switchBothSets(tmp);
}
// =======================================================================================
// Class ChromosomeProfileEndEntry
// =======================================================================================
// Stores the lastWindow and the insert size deviation of a specific insert. The ID equals the index of the 
// corresponding entry in the CyliclStartEntryTable.
struct ChromosomeProfileEndEntry
{
    TEntryIdType id;
    __uint32  lastWindow;

    ChromosomeProfileEndEntry() :
        id(0), lastWindow(0){}

    ChromosomeProfileEndEntry(TEntryIdType id, __int32 lastWindow) :
        id(id), lastWindow(lastWindow){}
};
// =======================================================================================
// Operator== overloads
// =======================================================================================
bool operator==(const ChromosomeProfileEndEntry & l, const ChromosomeProfileEndEntry & r)
{
    return (l.id == r.id && l.lastWindow == r.lastWindow);
}
bool operator==(const std::vector<ChromosomeProfileEndEntry> & l, const String<ChromosomeProfileEndEntry> & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
bool operator==(const String<std::vector<ChromosomeProfileEndEntry> > & l,
                const String<std::vector<ChromosomeProfileEndEntry> > & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
// =======================================================================================
// Operator< overload
// =======================================================================================
// Needed for std::upper_bound()
bool operator <(unsigned l, const ChromosomeProfileEndEntry & r)
{
    return l < r.lastWindow;
}
// =======================================================================================
// Class EndEntrySubset
// =======================================================================================
// Similar to class StartEntrySubset.
struct EndEntrySubset
{
    std::vector<ChromosomeProfileEndEntry>  subset;
    unsigned                                right;         // exclusiv right border.
    bool                                    readable;      // True to indicate that the subset is open for reading.

    EndEntrySubset():
    right(0), readable(false){}

    EndEntrySubset(unsigned r) :
    right(r), readable(false){}
    // =======================================================================================
    // Function clear()
    // =======================================================================================
    inline void clear()
    {
        subset.clear();
        right = 0;
        readable = false;
    }
    // =======================================================================================
    // Function reset()
    // =======================================================================================
    inline void reset(const unsigned & rightBorder)
    {
        subset.clear();
        right = rightBorder;
        readable = false;
    }
    // =======================================================================================
    // Function permitReading()
    // =======================================================================================
    inline void permitReading()
    {
        readable = true;
    }
    // =======================================================================================
    // Function forbidReading()
    // =======================================================================================
    inline void forbidReading()
    {
        readable = false;
    }
    inline ChromosomeProfileEndEntry at(const unsigned & i) const
    {
        SEQAN_ASSERT_LT(i, length(subset));
        return subset[i];
    }
};
inline unsigned length(const EndEntrySubset & e)
{
    return length(e.subset);
}
// =======================================================================================
// Operator== overloads
// =======================================================================================
bool operator==(const EndEntrySubset & l, const EndEntrySubset & r)
{
    if (length(l.subset) != length(r.subset))
        return false;
    if (l.right != r.right)
        return false;
    if (l.readable != r.readable)
        return false;
    Iterator<const std::vector<ChromosomeProfileEndEntry> >::Type  lIt  = begin(l.subset, Standard());
    Iterator<const std::vector<ChromosomeProfileEndEntry> >::Type lEnd = end(l.subset, Standard());
    Iterator<const std::vector<ChromosomeProfileEndEntry> >::Type rIt  = begin(r.subset, Standard());
    while (lIt != lEnd)
    {
        if (!(*lIt == *rIt))
            return false;
        ++lIt;
        ++rIt;
    }
    return true;
}
bool operator==(const String<EndEntrySubset> & l, const String<EndEntrySubset> & r)
{
    if (length(l) != length(r))
        return false;
    Iterator<const String<EndEntrySubset> >::Type lIt  = begin(l, Standard());
    Iterator<const String<EndEntrySubset> >::Type lEnd = end(l, Standard());
    Iterator<const String<EndEntrySubset> >::Type rIt  = begin(r, Standard());
    while (lIt != lEnd)
    {
        if (!(*lIt == *rIt))
            return false;
        ++lIt;
        ++rIt;
    }
    return true;
}
// =======================================================================================
// Overload for iterator function atBegin()
// =======================================================================================
inline bool atBegin(const Iterator<std::vector<ChromosomeProfileEndEntry>, Rooted>::Type & it)
{
    return it == begin(container(it), Rooted());
}
// =======================================================================================
// Class CyclicEndEntryTable
// =======================================================================================
// Similar to CyclicStartEntryTable. 
// A notable difference is that while the CylicStartEntryTable can automatically switch read and write sets,
// the CyclicEnd entry table cannot do so. Each time the sets in the CylicStartEntryTable are switched, 
// switchReadSets() has to be called for the CyclicEndEntryTable, then all procssing of the entries has to be done
// (including potential window merging) and then switchWriteSets has to be call for the CyclicEndEntryTable.
// After that, new elements can be added.
struct CyclicEndEntryTable
{
    String<EndEntrySubset>                                           sets;
    Iterator<std::vector<ChromosomeProfileEndEntry>, Rooted>::Type   nextRead;    // Element of next read action.
    unsigned                                                         dNext;       // Index in the current write set, where the count of active reads is reduced.
    unsigned                                                         dNextOffset; // Offset after look-ahead.
    unsigned                                                         writeSet;    // 0: first, 1: second, 2:third.
    unsigned                                                         endPosSet;   // EndPositions may lag behind write set.
    unsigned                                                         readSet;
    unsigned                                                         numWindows;  // Max number of windows per set.

    CyclicEndEntryTable(unsigned firstWin = 0, unsigned windows = 10000)
    {
        SEQAN_ASSERT_GT(windows, 0u);
        numWindows = windows;
        writeSet = 0;
        resize(sets, 3, Exact());
        sets[0].reset(firstWin + windows);
        sets[1].reset(sets[0].right + windows);
        sets[2].reset(sets[1].right + windows);
        readSet = 2;
        nextRead = begin(sets[readSet].subset, Rooted());
        dNext = 0;
        dNextOffset = 0;
        endPosSet = 0;
    }
    // =======================================================================================
    // Function needsSwitch()
    // =======================================================================================
    // Return true, if the given lastWindow requires to switch to a higher set.
    inline bool needsSwitch(const unsigned & lastWindow) const
    {
        if (lastWindow >= sets[writeSet].right)
            return true;
        else
            return false;
    }
    // =======================================================================================
    // Function switchWriteSet()
    // =======================================================================================
    // Switch between the three sets. Update the borders and the bool indicating the active set.
    // reset the set after the next set.
    inline void switchWriteSet()
    {
        EndEntrySubset & currentSet = sets[writeSet];
        EndEntrySubset & nextSet = sets[(writeSet + 1) % 3];
        EndEntrySubset & nextNextSet = sets[(writeSet + 2) % 3];
        nextNextSet.reset(nextSet.right + numWindows);  //We do not reset nextSet but the one after
        writeSet = (writeSet + 1) % 3;                  //Because nextSet might already contain entries.
        currentSet.permitReading(); // No neccessary if switchReadSets has been applied before.
    }
    inline void correctConsecutiveSwitch(unsigned & activeLoad)
    {
        if (endPosSet != (writeSet + 2) % 3) // There have been two consecutive switches.
        {
            endPosSet = writeSet;
            dNext = 0;
            SEQAN_ASSERT_EQ(dNextOffset, 0u);
            activeLoad = 0;
        }
        else if (length(sets[endPosSet]) == 0u)
        {
            SEQAN_ASSERT_EQ(dNext, 0u);
            SEQAN_ASSERT_EQ(activeLoad, 0u);
            endPosSet = writeSet;
        }
    }
    // =======================================================================================
    // Function switchReadSet()
    // =======================================================================================
    // Set nextRead to the begin of the next set.
    // Return false if the nextSet is empty, true otherwise.
    inline bool switchReadSet()
    {
        readSet = (readSet + 1) % 3;
        sets[readSet].permitReading();
        nextRead = begin(sets[readSet].subset, Rooted());
        if (atEnd(nextRead))
            return false;
        else
            return true;
    }
    // =======================================================================================
    // Function backwardSwitchReadSet()
    // =======================================================================================
    // Set nextRead to the end-1 of the next set.
    // Return false if the set is empty, true otherwise.
    inline bool backwardSwitchReadSet()
    {
        readSet = (readSet + 2) % 3;
        SEQAN_ASSERT(sets[readSet].readable);
        nextRead = end(sets[readSet].subset, Rooted());
        if (nextRead == begin(container(nextRead), Rooted()))
        {
            return false;
        }
        else
        {
            --nextRead;
            return true;
        }
    }
    // =======================================================================================
    // Function switchBothSets()
    // =======================================================================================
    // Calls f() and switchReadSet().
    // This should alwasy result in nextRead being set to the beginning of the previous writeSet.
    // Eg. readSet == 1, writeSet == 2  -> readSet == 2; writeSet == 0;
    // Return false if the next set for reading is empty, true otherwise.
    inline bool switchBothSets()
    {
        switchWriteSet();
        return switchReadSet();
    }
    // =======================================================================================
    // Function add()
    // =======================================================================================
    // Add the entry to the correspondig subset.
    inline void add(const unsigned & id, const unsigned & lastWin)
    {
        unsigned whichSet = writeSet;
        if (needsSwitch(lastWin))
        {
            whichSet = (whichSet + 1) % 3;
        }
        SEQAN_ASSERT_EQ(sets[whichSet].readable, false);    // Never insert into a readable subset.
        std::vector<ChromosomeProfileEndEntry>& currentSubset = sets[whichSet].subset;
        std::vector<ChromosomeProfileEndEntry>::iterator at = std::upper_bound(currentSubset.begin(),
                                                                               currentSubset.end(),
                                                                               lastWin);
        currentSubset.insert(at, ChromosomeProfileEndEntry(id, lastWin));
    }
    // =======================================================================================
    // Function goNext()
    // =======================================================================================
    // Move the iterator for reading forward by one.
    // Return false if the the iterator reached the end of the current subset and requires a switch.
    inline bool goNext()
    {
        SEQAN_ASSERT(!atEnd(nextRead));
        ++nextRead;
        if (atEnd(nextRead))
            return false;
        else
            return true;
    }
    // =======================================================================================
    // Function goPrevious()
    // =======================================================================================
    // Move the iterator for reading back by one.
    // Return true on succes and false if the set has been switch but the new set is empty.
    inline bool goPrevious()
    {
        SEQAN_ASSERT(sets[readSet].readable);
        if (atBegin(nextRead))
        {
            return backwardSwitchReadSet();
        }
        else
        {
            --nextRead;
            return true;
        }
    }
    // =======================================================================================
    // Function peekPreviousLastWindow()
    // =======================================================================================
    // Return the lastWindow of the previous position without changing the iterator.
    // Return 0 if the previous position is invalid.
    inline unsigned peekPreviousLastWindow() const
    {
        if (atBegin(nextRead))
        {
            unsigned prevSetNum = (readSet + 2) % 3;
            if (empty(sets[prevSetNum].subset) || sets[readSet].right < sets[prevSetNum].right)
                return 0;
            else
                return  back(sets[prevSetNum].subset).lastWindow;
        }
        else
        {
            return (nextRead - 1)->lastWindow;   
        }
    }
    // =======================================================================================
    // Function getEntry()
    // =======================================================================================
    // Return a reference to the element nextReas is pointing to.
    inline ChromosomeProfileEndEntry & getEntry() const
    {
        SEQAN_ASSERT(!atEnd(nextRead));
        return *nextRead;
    }
    // =======================================================================================
    // Function getEndCount()
    // =======================================================================================
    // Return the number of reads that end in the current writing set between the position of [dNext and pos[.
    // Also updates the position of dNext.
    inline unsigned getEndCount(const unsigned & pos)
    {
        bool lookAhead = false;                             // Indicates that we have to look in endPosSet + 1.
        unsigned i = dNext;
        unsigned c = 0;                                     // Counter for removed reads.
        unsigned cEndPosSet = endPosSet;
        unsigned n = length(sets[cEndPosSet]);
        if (i == n)
        {
            if (endPosSet != writeSet)  // We are done with the current endSet. Update it to the current writeSet.
            {
                endPosSet = (endPosSet + 1) % 3;
                dNext = dNextOffset;
                dNextOffset = 0;
                i = dNext;
                cEndPosSet = endPosSet;
                n = length(sets[cEndPosSet]);

            }
            else
            {
                return 0;
            }
        }
        while(sets[cEndPosSet].at(i).lastWindow < pos)
        {
            if (lookAhead)
                ++dNextOffset;
            ++c;
            ++i;
            if (i == n)
            {
                if (endPosSet != writeSet)  // We are done with the current endSet. Update it to the current writeSet.
                {
                    endPosSet = (endPosSet + 1) % 3;
                    dNext = dNextOffset;
                    dNextOffset = 0;
                    return c;
                }
                else
                {
                    // some of the previous end pos. have been added to the next writeSet,
                    // but the current endSet might still be used later.
                    // So we cannot yet permanently switch to the next endPosSet.
                    cEndPosSet = (cEndPosSet + 1) % 3;
                    n = length(sets[cEndPosSet]);
                    if (n == 0)
                        break;
                    i = dNextOffset;
                    lookAhead = true;
                }
            }
        }
        dNext += c - dNextOffset;
        return c;
    }
    // =======================================================================================
    // Function getId()
    // =======================================================================================
    // Return the ID of the element nextRead is pointing to.
    inline unsigned getId() const
    {
        SEQAN_ASSERT(!atEnd(nextRead));
        return nextRead->id;
    }
    // =======================================================================================
    // Function getLastWin()
    // =======================================================================================
    // Return lastWin of the element nextRead is pointing to.
    inline unsigned getLastWin() const
    {
        SEQAN_ASSERT(sets[readSet].readable);
        SEQAN_ASSERT(!atEnd(nextRead));
        return nextRead->lastWindow;
    }
};
// =======================================================================================
// Operator== overloads
// =======================================================================================
bool operator==(const CyclicEndEntryTable & l, const CyclicEndEntryTable & r)
{

    if (!(l.sets == r.sets))
        return false;
    if (position(l.nextRead) != position(r.nextRead))
        return false;
    if (l.readSet != r.readSet)
        return false;
    if (l.writeSet != r.writeSet)
        return false;
    if (l.numWindows != r.numWindows)
        return false;
    else
        return true;
}
bool operator==(const String<CyclicEndEntryTable> & l, const String<CyclicEndEntryTable> & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
// =======================================================================================
// Function switchReadSets()
// =======================================================================================
// Wrapper for updating the read permissions and calling switchReadSet on a String of CyclicEndEntryTable.
inline void switchReadSets(String<CyclicEndEntryTable> & s)
{
    for (Iterator<String<CyclicEndEntryTable> >::Type it = begin(s); it != end(s); ++it)
    {
        it->sets[(it->readSet + 1) % 3].permitReading();;
        it->switchReadSet();
    }
}
// =======================================================================================
// Function switchWriteSets()
// =======================================================================================
// Wrapper for calling switchWriteSet on a String of CyclicEndEntryTable.
inline void switchWriteSets(String<CyclicEndEntryTable> & s)
{
    for (Iterator<String<CyclicEndEntryTable> >::Type it = begin(s); it != end(s); ++it)
    {
        it->switchWriteSet();
    }
}
// =======================================================================================
// Function atBegin()
// =======================================================================================
// Return true if nextRead points to the beginning of its container.
inline bool atBegin(const CyclicEndEntryTable & c)
{
    return c.nextRead == begin(container(c.nextRead), Rooted());
}
// =======================================================================================
// Function atEnd()
// =======================================================================================
// Return true if nextRead points to the end its container.
inline bool atEnd(const CyclicEndEntryTable & c)
{
    return atEnd(c.nextRead);
}
// =======================================================================================
// Function switchSets()
// =======================================================================================
// Wrapper class for applying switchSet() on every CyclicEndEntryTable in the string.
// Only used for convenience during unit tests.
inline void switchSets(String<CyclicEndEntryTable> & s)
{
    for (Iterator<String<CyclicEndEntryTable> >::Type it = begin(s); it != end(s); ++it)
    {
        it->switchBothSets();
    }
}
// =======================================================================================
// Class ChromosomeProfile
// =======================================================================================
// Class holding and managing the CylcicStartEntryTable and CyclicEndEntryTable for each readGroup.
// Also maintains the sets of active reads, which contains the IDs to the entries in the CylcicStartEntryTable.
struct ChromosomeProfile
{
    typedef std::unordered_set<TEntryIdType>                                 TActiveSet;
    typedef Iterator<String<ChromosomeProfileStartEntry>, Rooted>::Type      TChromosomeProfileStartEntryIterator;
    typedef Iterator<std::vector<ChromosomeProfileEndEntry>, Rooted>::Type   TChromosomeProfileEndEntryIterator;

    String<CyclicStartEntryTable>                   startProfiles;  // One per read group
    String<CyclicEndEntryTable>                     endProfiles;    // One per read group
    String<TActiveSet>                              activeReads;    // Set of IDs of active reads, one per RG
    // TODO: Maybe add string containing the maximum of each window to speed up checkActiveReads().
    unsigned                                        globalMinPos;   // First position of the profile.
    unsigned                                        currentPos;     // Current position on chromosome (0-based);
    String<double>                          avgNewReadsPerWindow; // avg. # of new reads per window for each read Group.
    double                                  totalAvgNewReadsPerWindow; // sum of avgNewReadsPerWindow;
    bool                                    profilesAtEnd;      // That the start and endProfiles have ended.
    unsigned                                currentRightBorder; // Right border of the currently active profile subset.
                                                                // Corresponds to startProfiles[rg].sets[readSet].right.
    String<unsigned>                        activeLoad;// Number of active reads in current writing window. One per RG.
    String<unsigned>                        maxLoad; // while the load is above this threshold no new reads are added.

    ChromosomeProfile(unsigned numReadGroups, String<unsigned> & maximumLoad, unsigned bufferSize = 100)
    {
        resize(startProfiles, numReadGroups, CyclicStartEntryTable(0, bufferSize), Exact());
        resize(endProfiles, numReadGroups, CyclicEndEntryTable(0, bufferSize), Exact());
        resize(activeReads, numReadGroups, Exact());
        resize(avgNewReadsPerWindow, numReadGroups, 0.0, Exact());
        resize(activeLoad, numReadGroups, 0, Exact());
        globalMinPos = maxValue<unsigned>();
        currentPos = maxValue<unsigned>();
        profilesAtEnd = false;
        currentRightBorder = startProfiles[0].sets[2].right;
        move(maxLoad, maximumLoad);
    }
    // =======================================================================================
    // Function add()
    // =======================================================================================
    // Append values to the read group's profile. Exspects the added Elements to be ordered by firstWindow.
    // If necessary, one switch is performed, but not more.
    // Return 0 if the element has been added without a switch occuring in the startProfile.
    // Return 1 if the element has been added and a switch occured in the startProfile.
    // Return 2 if the element has not been added and a switch occured in the startProfile.
    inline void add(unsigned readGroup, unsigned firstWindow, unsigned lastWindow, int deviation)
    {
        SEQAN_ASSERT_LEQ(readGroup, length(startProfiles));
        //The starting positions are already sorted by startingWindow.
        CyclicStartEntryTable& currentStartEntryTable = startProfiles[readGroup];
        CyclicEndEntryTable&   currentEndEntryTable = endProfiles[readGroup];
        unsigned id = currentStartEntryTable.insertionCount;
        if (activeLoad[readGroup] >= maxLoad[readGroup])
        {
            // Don't add new reads, only look for closing reads.
            unsigned closingReads = currentEndEntryTable.getEndCount(firstWindow);
            SEQAN_ASSERT_GT(activeLoad[readGroup], closingReads);
            activeLoad[readGroup] -= closingReads;
            if (activeLoad[readGroup] < maxLoad[readGroup]) // Maybe it dropped below the threshold now.
            {
                currentStartEntryTable.add(firstWindow, deviation);
                currentEndEntryTable.add(id, lastWindow);
                ++activeLoad[readGroup];
            }
        }
        else
        {
            currentStartEntryTable.add(firstWindow, deviation);
            currentEndEntryTable.add(id, lastWindow);
            ++activeLoad[readGroup];
            unsigned closingReads = currentEndEntryTable.getEndCount(firstWindow);
            SEQAN_ASSERT_GT(activeLoad[readGroup], closingReads);
            activeLoad[readGroup] -= closingReads;
        }
       // std::cout << "RG:" << readGroup << "\tPos:" << firstWindow << "\tLoad:" << activeLoad[readGroup] << std::endl;
    }
    // =======================================================================================
    // Function initializeActiveReads()
    // =======================================================================================
    // Call initializeIterators() and add the first reads to the set of active reads
    // and advance iterators for startWindows to the next entry where appropriate (not for endWindows!).
    // This functions should only be called after the first add() returned a value > 0.
    inline void initializeActiveReads()
    {
        for (__uint32 rg = 0; rg < length(startProfiles); ++rg)        // Find the minimum position.
        {
            CyclicStartEntryTable& currentStartEntryTable = startProfiles[rg];
            currentStartEntryTable.nextRead = begin(currentStartEntryTable.sets[currentStartEntryTable.readSet].subset, Rooted());
            if (atEnd(startProfiles[rg]))                              // Skip read group, in case it is empty.
                continue;
            unsigned pos = startProfiles[rg].getFirstWin();
            if (pos < currentPos)
            {
                currentPos = pos;
            }
        }
        globalMinPos = currentPos;
        Iterator<String<CyclicStartEntryTable>, Rooted>::Type currentProfileIt(begin(startProfiles));
        Iterator<String<TActiveSet>, Standard >::Type itA(begin(activeReads, Standard()));
        while(!atEnd(currentProfileIt))
        {
            CyclicStartEntryTable& currentStartEntryTable = *currentProfileIt;
            if (!atEnd(currentStartEntryTable) && itA->empty())
            {
                unsigned pos = currentStartEntryTable.getFirstWin();
                if (pos != currentPos)      // Skip this read group, if its first window is != the global first window.
                {
                    ++currentProfileIt;
                    ++itA;
                    continue;
                }
                TEntryIdType id = currentProfileIt->sets[currentProfileIt->readSet].idOffset;
                (*itA).emplace(id);
                ++id;
                bool endOfSubSet = !currentStartEntryTable.goNext();
                while(!endOfSubSet && currentStartEntryTable.getFirstWin() == pos)
                {
                    (*itA).emplace(id);
                    ++id;
                    endOfSubSet = !currentStartEntryTable.goNext();
                }
            }
            ++currentProfileIt;
            ++itA;
        }
    }
    // =======================================================================================
    // Function nextStartWindow()
    // =======================================================================================
    // Advance to next window in the RG's string of start entries and add IDs to its active set.
    // Return true on success and false, if the last window of the subset has been processed.
    inline unsigned nextStartWindow(const __uint32 & rg)
    {
        CyclicStartEntryTable& currentStartEntryTable =  startProfiles[rg];
        SEQAN_ASSERT(!atEnd(currentStartEntryTable));
        TActiveSet& currentActiveSet = activeReads[rg];
        unsigned currentWindow = currentStartEntryTable.getFirstWin();
        bool noSwitchNecessary = true;
        while (noSwitchNecessary && currentWindow == currentStartEntryTable.getFirstWin())
        {
            currentActiveSet.emplace(currentStartEntryTable.getEntryId());
            noSwitchNecessary = currentStartEntryTable.goNext();
        }
        return noSwitchNecessary;
    }
    // =======================================================================================
    // Function nextEndWindow()
    // =======================================================================================
    // Advance to next window in the RG's string of end entries and remove IDs from its active set.
    // Return true on success and false, if the last window of a set has been processed.
    inline bool nextEndWindow(const __uint32 & rg)
    {
        CyclicEndEntryTable& currentEndEntryTable = endProfiles[rg];
        if (atEnd(currentEndEntryTable))   //TODO: This is strange... check
        {
            currentEndEntryTable.goNext();
            if (atEnd(currentEndEntryTable))    // Should never happen. TODO: test and remove.
                return false;
        }
        TActiveSet& currentActiveSet = activeReads[rg];
        unsigned currentWindow = currentEndEntryTable.getLastWin();
        while(currentWindow == currentEndEntryTable.getLastWin())
        {
            currentActiveSet.erase(currentEndEntryTable.getId());
            if (!currentEndEntryTable.goNext())
                return false;
        }
        return true;
    }
    // =======================================================================================
    // Function nextWindow()
    // =======================================================================================
    // Take the windowShift parameter and call nextStartWindow() and nextEndWindow for each read group
    //  if appropriate, thus updating the iterators to the start- and end-entries and updating the set of active reads.
    // Return false the subsets of ALL readGroups require a switch and no shift has been performed.
    inline bool nextWindow(const unsigned & shift)
    {
        SEQAN_ASSERT_GT(shift, 0u);
        bool someSubsetStillGood = false;
        currentPos += shift;
        if (currentPos >= currentRightBorder)
        {
            currentPos -= shift;
            return false;
        }
        for (__uint32 rg = 0; rg < length(startProfiles); ++rg)
        {
            CyclicStartEntryTable& currentStartEntryTable =  startProfiles[rg];
            CyclicEndEntryTable& currentEndEntryTable =  endProfiles[rg];
            if (!(atEnd(currentEndEntryTable) && atEnd(currentStartEntryTable)))
                someSubsetStillGood = true;
            while (!atEnd(currentStartEntryTable))
            {
                if (currentStartEntryTable.getFirstWin() <= currentPos)
                    nextStartWindow(rg);
                else
                    break;
            }
            while (!atEnd(currentEndEntryTable))
            {
                if (currentEndEntryTable.getLastWin() < currentPos) // < because the last window is still valid.
                    nextEndWindow(rg);
                else
                    break;
            }
        }
        if (!someSubsetStillGood)
            currentPos -= shift;
        return (someSubsetStillGood);
    }
    // =======================================================================================
    // Function previousStartWindow()
    // =======================================================================================
    // Move back one window in the RG's string of start entries and remove IDs from its active set.
    // Basically undos all changes made by one call of nextStartWindow.
    // Return true on success and false, if the first window has been reached.
    inline bool previousStartWindow(const __uint32 & rg)
    {
        CyclicStartEntryTable& currentStartEntryTable =  startProfiles[rg];
        TActiveSet& currentActiveSet = activeReads[rg];
        unsigned currentWindow = currentStartEntryTable.peekPreviousFirstWindow();
        if (currentWindow <= globalMinPos)
            return false;
        bool goOn;
        unsigned previousWindow;
        do
        {
            goOn = currentStartEntryTable.goPrevious();
            currentActiveSet.erase(currentStartEntryTable.getEntryId());
            previousWindow = currentStartEntryTable.peekPreviousFirstWindow();
        }
        while(goOn && (currentWindow == previousWindow) &&  (previousWindow > globalMinPos));
        return true;
    }
    // =======================================================================================
    // Function previousEndWindow()
    // =======================================================================================
    // Move back one window in the RG's string of end entries and add IDs to its active set.
    // Basically undos all changes made by one call of nextEndWindow.
    // Return false if the set has been switched, true otherwise.
    inline bool previousEndWindow(const __uint32 & rg)
    {
        CyclicEndEntryTable& currentEndEntryTable =  endProfiles[rg];
        TActiveSet& currentActiveSet = activeReads[rg];
        bool sameSet = currentEndEntryTable.goPrevious();
        if(empty(container(currentEndEntryTable.nextRead)))
            return false;
        unsigned currentWindow = currentEndEntryTable.getLastWin();
        while(currentWindow == currentEndEntryTable.getLastWin())
        {
            currentActiveSet.emplace(currentEndEntryTable.getId());
            unsigned prevWindow = currentEndEntryTable.peekPreviousLastWindow();
            if (prevWindow  == currentWindow)
                currentEndEntryTable.goPrevious();
            else
                return sameSet;
        }
        return sameSet;
    }
    // =======================================================================================
    // Function previousWindow()
    // =======================================================================================
    // Take the windowShift parameter and call previousStartWindow() and previousEndWindow for each read group
    //  if appropriate, thus undoing calls of nextWindow.
    // Return true if at least one shift could be performed, false otherwise
    inline bool previousWindow(const unsigned & shift)
    {
        if (currentPos - shift < globalMinPos)
            return false;
        currentPos -= shift;
        bool shifted = false;
        for (__uint32 rg = 0; rg < length(startProfiles); ++rg)
        {
            CyclicEndEntryTable& currentEndEntryTable =  endProfiles[rg];
            unsigned currentSetStart = currentEndEntryTable.sets[currentEndEntryTable.readSet].right -
                                       currentEndEntryTable.numWindows;
            if (empty(container(currentEndEntryTable.nextRead)))
            {
                if (currentSetStart >= currentPos)
                    shifted &= previousEndWindow(rg);
            }
            else if (atBegin(currentEndEntryTable))
            {
                unsigned prevSetBorder = currentEndEntryTable.sets[(currentEndEntryTable.readSet + 2) % 3].right;
                if (currentEndEntryTable.getLastWin() >= currentPos && prevSetBorder == currentSetStart)
                    shifted &= previousEndWindow(rg);
            }
            else if (currentEndEntryTable.peekPreviousLastWindow() >= currentPos)
            {
                shifted &= previousEndWindow(rg);
            }
            CyclicStartEntryTable& currentStartEntryTable =  startProfiles[rg];
            if (empty(container(currentStartEntryTable.nextRead)))
            {
                if (currentSetStart >= currentPos)
                    shifted &= previousStartWindow(rg);
            }
            else if (atBegin(currentStartEntryTable))
            {
                if (currentStartEntryTable.getFirstWin() > currentPos)
                    shifted &= previousStartWindow(rg);
            }
            else if (currentStartEntryTable.peekPreviousFirstWindow() > currentPos)
            {
                    shifted &= previousStartWindow(rg);
            }
        }
        return true;
    }
    // =======================================================================================
    // Function mergeWindows()
    // =======================================================================================
    // Merges the next 'numWindows' windows for all read groups.
    // This is done by repeating the standard update process of nextWindow, without updating the endIteraros,
    // thus not removing any entries from the active set.
    // This function changes the read and write set entries for startProfiles and endProfiles, leaving them in 
    // an incosistent state (endProfile being left behing). Befor continuing, this has to be fixed by calling
    // endIteratorCatchUp().
    // Return true if at least one shift could be performed, false otherwise
    inline bool mergeWindows(const unsigned & numWindows, const unsigned & shift)
    {
        if (numWindows < 2)
            return false;
        bool end = true;
        currentPos += shift * (numWindows - 1);
        for (__uint32 rg = 0; rg < length(startProfiles); ++rg)
        {
            CyclicStartEntryTable& currentStartEntryTable =  startProfiles[rg];
            for (unsigned i = 1; i < numWindows; ++i)
            {
                if (atEnd(currentStartEntryTable) &&
                    currentStartEntryTable.sets[(currentStartEntryTable.readSet + 1) % 3].readable)
                {
                    currentStartEntryTable.switchReadSet();
                }
                if (!atEnd(currentStartEntryTable))
                {
                    end = false;
                    if (currentStartEntryTable.getFirstWin() <= currentPos)
                    {
                        if (!nextStartWindow(rg))
                            continue;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        return (!end);
    }
    // =======================================================================================
    // Function endIteratorCatchUp()
    // =======================================================================================
    // After mergeWindows() has been called and the window has been processed, the active set has to be brought
    //  into an state, as if the merge never had happend. This is done by catching up on the missing calls of
    //  next endWindow()
    // Return true if at least one update could be performed, false otherwise.
    inline bool endIteratorCatchUp(const unsigned & numWindows)
    {
        if (numWindows < 2)
            return false;
        bool end = true;
        for (__uint32 rg = 0; rg < length(endProfiles); ++rg)
        {
            CyclicEndEntryTable& currentEndEntryTable =  endProfiles[rg];
            for (unsigned i = 1; i < numWindows; ++i)
            {
                if (atEnd(currentEndEntryTable) &&
                    currentEndEntryTable.sets[(currentEndEntryTable.readSet + 1) % 3].readable)
                {
                    currentEndEntryTable.switchReadSet();
                }
                if (!atEnd(currentEndEntryTable))
                {
                    end = false;
                    if (currentEndEntryTable.getLastWin() < currentPos)         // < because last window is still valid.
                    {
                        if (!nextEndWindow(rg))
                            continue;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        return (!end);
    }
    // =======================================================================================
    // Function getActiveReadsDeviations()
    // =======================================================================================
    // Get the insert size deviations of the currently active reads for all given read groups.
    // Return the number of active reads in the given read groups.
    inline unsigned getActiveReadsDeviations(String<int> & deviations,
                                             const TReadGroupIndices & rgs,
                                             unsigned minCov) const
    {
        SEQAN_ASSERT(empty(deviations));
        unsigned readCount = 0;
        Iterator<const TReadGroupIndices >::Type cReadGroupsIt(begin(rgs));
        Iterator<const TReadGroupIndices >::Type cReadGroupsItEnd(end(rgs));
        while (cReadGroupsIt != cReadGroupsItEnd)               // For each active read group...
        {
            readCount += activeReads[*cReadGroupsIt].size();    //...get the total number of active reads.
            ++cReadGroupsIt;
        }
        if (readCount < minCov)
            return readCount;
        resize(deviations, readCount);
        unsigned i = 0;
        for (cReadGroupsIt = begin(rgs); cReadGroupsIt != cReadGroupsItEnd; ++cReadGroupsIt)
        {
            TActiveSet::const_iterator where(activeReads[*cReadGroupsIt].begin());
            TActiveSet::const_iterator whereEnd(activeReads[*cReadGroupsIt].end());
            while (where != whereEnd)
            {
                deviations[i] = startProfiles[*cReadGroupsIt].getDeviationAt(*where);
                ++i;
                ++where;
            }
        }
        return readCount;
    }
    // =======================================================================================
    // Function getActiveReadsSizes()
    // =======================================================================================
    // Get the total insert sizes of the currently active reads for all given read groups.
    // Works like getActiveReadsDeviations, but adds the median of the read group to the deviations.
    inline unsigned getActiveReadsSizes(String<int> & insertSizes,
                                        const String<Histogram> & hists,
                                        const TReadGroupIndices & rgs,
                                        unsigned minCov = 1) const
    {
        SEQAN_ASSERT(empty(insertSizes));
        unsigned readCount = 0;
        Iterator<const TReadGroupIndices >::Type cReadGroupsIt(begin(rgs));
        Iterator<const TReadGroupIndices >::Type cReadGroupsItEnd(end(rgs));
        while (cReadGroupsIt != cReadGroupsItEnd)               // For each active read group...
        {
            readCount += activeReads[*cReadGroupsIt].size();    //...get the total number of active reads.
            ++cReadGroupsIt;
        }
        if (readCount < minCov)
            return readCount;
        resize(insertSizes, readCount);
        unsigned i = 0;
        for (cReadGroupsIt = begin(rgs); cReadGroupsIt != cReadGroupsItEnd; ++cReadGroupsIt)
        {
            unsigned currentMedian = hists[*cReadGroupsIt].median;
            TActiveSet::const_iterator where(activeReads[*cReadGroupsIt].begin());
            TActiveSet::const_iterator whereEnd(activeReads[*cReadGroupsIt].end());
            while (where != whereEnd)
            {
                insertSizes[i] = startProfiles[*cReadGroupsIt].getDeviationAt(*where) + currentMedian;
                ++i;
                ++where;
            }
        }
        return readCount;
    }
    // =======================================================================================
    // Function getActiveReadsFirstLast()
    // =======================================================================================
    // Return the lowest firstWin and highest lastWin of the currently active reads for the given read groups.
    inline Pair<unsigned> getActiveReadsFirstLast(const String<Histogram> & hists, const TReadGroupIndices & rgs) const
    {
        Iterator<const TReadGroupIndices >::Type cReadGroupsIt(begin(rgs));
        Iterator<const TReadGroupIndices >::Type cReadGroupsItEnd(end(rgs));
        unsigned currentMin = maxValue<unsigned>();
        unsigned currentMax = 0;
        for (cReadGroupsIt = begin(rgs); cReadGroupsIt != cReadGroupsItEnd; ++cReadGroupsIt)
        {
            const Histogram & hist = hists[*cReadGroupsIt];
            int currentMedian = hist.median;
            int currentDoubleReadLength = 2 * hist.readLength;
            for (TActiveSet::const_iterator where(activeReads[*cReadGroupsIt].begin());
                 where != activeReads[*cReadGroupsIt].end();
                 ++where)
            {
                const unsigned firstWin = startProfiles[*cReadGroupsIt].getFirstWinAt(*where);
                int insertSize = startProfiles[*cReadGroupsIt].getDeviationAt(*where) + currentMedian;
                insertSize = std::max(0, insertSize - currentDoubleReadLength);
                const unsigned lastWin = ((firstWin + insertSize) / 30) * 30; // TODO: As parameter?
                if (firstWin < currentMin)
                    currentMin = firstWin;
                if (lastWin > currentMax)
                    currentMax = lastWin;
            }
        }
        if (currentMin == maxValue<unsigned>())
            currentMin = 0;
        return Pair<unsigned>(currentMin, currentMax);
    }
    // =======================================================================================
    // Function updateSupportFirstLast()
    // =======================================================================================
    // Similar to getActiveReadsFirstLast(), but returns the values only for reads between lower and upper (inlcusive).
    // Checks for the lowest first and the highest last window.
    // Does not return the value, but updates the given pair.
    // TODO: This method might be sensitive to samples with bad mappings. Make it more robust by not simply taking
    //   the min but some kind of weighted value or a mayority vote.
    inline void updateSupportFirstLast(Pair<unsigned> & suppFirstLast,
                                       const String<Histogram> & hists,
                                       const TReadGroupIndices & rgs,
                                       const int & lower,
                                       const int & upper) const
    {
        Iterator<const TReadGroupIndices >::Type cReadGroupsIt(begin(rgs));
        Iterator<const TReadGroupIndices >::Type cReadGroupsItEnd(end(rgs));
        for (cReadGroupsIt = begin(rgs); cReadGroupsIt != cReadGroupsItEnd; ++cReadGroupsIt)
        {
            const Histogram & hist = hists[*cReadGroupsIt];
            int currentMedian = hist.median;
            int currentDoubleReadLength = 2 * hist.readLength;
            for (TActiveSet::const_iterator where(activeReads[*cReadGroupsIt].begin());
                 where != activeReads[*cReadGroupsIt].end();
                 ++where)
            {
                int deviation = startProfiles[*cReadGroupsIt].getDeviationAt(*where);
                if (deviation >= lower && deviation <= upper)
                {
                    int insertSize = std::max(0, deviation + currentMedian - currentDoubleReadLength);
                    const unsigned firstWin = startProfiles[*cReadGroupsIt].getFirstWinAt(*where);
                    const unsigned lastWin = ((firstWin + insertSize) / 30) * 30; // TODO: As parameter?
                    if (firstWin > suppFirstLast.i1)
                        suppFirstLast.i1 = firstWin;
                    if (lastWin < suppFirstLast.i2)
                        suppFirstLast.i2 = lastWin;
                }
            }
        }
    }
    // =======================================================================================
    // Function getActiveReadsNum()
    // =======================================================================================
    // Return the number or currently active read pairs for the given rg.
    inline unsigned getActiveReadsNum(const TEntryIdType & rg)
    {
        return (activeReads[rg].size());
    }
    //calls getActiveReadsNum() for a string of read groups.
    inline unsigned getActiveReadsNum(const TReadGroupIndices & rgs)
    {
        unsigned  c = 0;
        for (unsigned i = 0; i < length(rgs); ++i)
            c += getActiveReadsNum(rgs[i]);
        return c;
    }
    // Return the number or currently starting active read pairs for the given rg with insert-size deviation >= minlen.
    inline unsigned getActiveReadsNum(const TEntryIdType & rg, const unsigned & currentWindow, const int & minLen)
    {
        TActiveSet::const_iterator it(activeReads[rg].begin());
        TActiveSet::const_iterator itEnd(activeReads[rg].end());
        unsigned  c = 0;
        while(it != itEnd)
        {
            const ChromosomeProfileStartEntry & currentEntry = startProfiles[rg][*it];
            if (currentEntry.insertSizeDeviation >= minLen && currentEntry.firstWindow == currentWindow)
                ++c;
            ++it;
        }
        return c;
    }
    // Return the number or currently starting active read pairs for the all rg's with insert-size >= minlen.
    inline unsigned getActiveReadsNum(const TReadGroupIndices & rgs,
                                      const unsigned & currentWindow,
                                      const int & minLen)
    {
        unsigned  c = 0;
        for (unsigned i = 0; i < length(rgs); ++i)
            c += getActiveReadsNum(rgs[i], currentWindow, minLen);
        return c;
    }
    // =======================================================================================
    // Function checkActiveReads()
    // =======================================================================================
    // Return true if the current activeSet contains read pairs whose deviation is within i1 and i2 of borders.
    // Return false otherwise.
    inline bool checkActiveReads(const TEntryIdType & rg, const Pair<__int32> & borders) const
    {
        TActiveSet::const_iterator it(activeReads[rg].begin());
        TActiveSet::const_iterator itEnd(activeReads[rg].end());
        while(it != itEnd)
        {
            const ChromosomeProfileStartEntry & currentEntry = startProfiles[rg].getEntryAt(*it);
            if (currentEntry.insertSizeDeviation >= borders.i1 && currentEntry.insertSizeDeviation <= borders.i2)
                return true;
            ++it;
        }
        return false;
    }
    // Overload of checkActiveReads() for application on muliple read groups.
    // Return true, if checkActiveReads reads is true for at least one RG.
    inline bool checkActiveReads(const TRGs & rgs, const String<Pair<__int32> >& borders) const
    {
        Iterator<const TRGs, Rooted>::Type sampleIt = begin(rgs, Rooted());
        while (!atEnd(sampleIt))
        {
            Iterator<const TReadGroupIndices, Rooted>::Type rgIt = begin(*sampleIt, Rooted());
            while(!atEnd(rgIt))
            {
                if (checkActiveReads(*rgIt, borders[*rgIt]))
                    return true;
                ++rgIt;
            }
            ++sampleIt;
        }
        return false;
    }
    // =======================================================================================
    // Function getSingleDeviation()
    // =======================================================================================
    inline int getSingleDeviation(const TEntryIdType & rg,
                                  const TActiveSet::const_iterator & actReadsIt) const
    {
        return startProfiles[rg].getDeviationAt(*actReadsIt);
    }
    // =======================================================================================
    // Function goToPosition()
    // =======================================================================================
    // Call nextWindow or previousWindow with 'shift' until 'tagetPos' has been reached
    //  or - in case of an incompatible 'shift' - passed
    inline void goToPosition(const unsigned & targetPos, const unsigned & shift)
    {
        if (currentPos > targetPos)
        {
            while(currentPos > targetPos)
            {
                previousWindow(shift);
            }
        }
        else
        {
            while(currentPos < targetPos)
            {
                nextWindow(shift);
            }
        }
    }
    // =======================================================================================
    // Function resetTo()
    // =======================================================================================
    // Reset the start and end profiles, s.t. the right border of all first sets is pos + numwindows.
    inline void resetTo(unsigned pos)
    {
        unsigned r1 = pos + startProfiles[0].numWindows;
        unsigned r2 = pos + 2 * startProfiles[0].numWindows;
        unsigned r3 = pos + 3 * startProfiles[0].numWindows;
        for (Iterator<String<CyclicStartEntryTable> >::Type it = begin(startProfiles); it != end(startProfiles); ++it)
        {
            unsigned idO = it->insertionCount;
            it->sets[0].reset(idO, r1);
            it->sets[1].reset(idO, r2);
            it->sets[2].reset(idO, r3);
            it->writeSet = 0;
            it->readSet = 2;
            it->nextRead = begin(it->sets[2].subset);
        }
        for (Iterator<String<CyclicEndEntryTable> >::Type it = begin(endProfiles); it != end(endProfiles); ++it)
        {
            it->sets[0].reset(r1);
            it->sets[1].reset(r2);
            it->sets[2].reset(r3);
            it->writeSet = 0;
            it->readSet = 2;
            it->nextRead = begin(it->sets[2].subset, Rooted());
        }
        currentRightBorder = startProfiles[0].sets[2].right;
    }
    // =======================================================================================
    // Function fullReset()
    // =======================================================================================
    // Completely reset the start and end profiles.
    inline void fullReset()
    {
        for (Iterator<String<CyclicStartEntryTable> >::Type it = begin(startProfiles); it != end(startProfiles); ++it)
           it->insertionCount = 0;
        for (Iterator<String<double> >::Type it = begin(avgNewReadsPerWindow); it != end(avgNewReadsPerWindow); ++it)
            *it = 0.0;
        resetTo(0);
        globalMinPos = maxValue<unsigned>();
        currentPos = maxValue<unsigned>();
        profilesAtEnd = false;
        currentRightBorder = startProfiles[0].sets[2].right;
    }
    inline void reportStatus()
    {
        std::cout << "ChromosomeProfile status report" << std::endl;
        std::cout << "currentPos:\t\t"  << currentPos
                  << "\nglobalMinPos:\t\t" << globalMinPos
                  << "\nprofilesAtEnd:\t\t" << profilesAtEnd
                  << "\ncurrentRightBorder:\t" << currentRightBorder
                  << "\n#RG's:\t\t\t" << length(activeReads)
                  << "\n#StartProf.:\t\t" << length(startProfiles)
                  << "\n#EndProf.:\t\t" << length(endProfiles)
                  << "\nStart read set:\t\t" << startProfiles[0].readSet
                  << "\nStart write set:\t" << startProfiles[0].writeSet
                  << "\nStart insert #:\t\t" << startProfiles[0].insertionCount
                  << "\nEnd read set:\t\t"  << endProfiles[0].readSet
                  << "\nEnd write set:\t\t" << endProfiles[0].writeSet << std::endl;
    }
};
// =======================================================================================
// Operator== overloads
// =======================================================================================
inline bool operator==(const String<ChromosomeProfile::TChromosomeProfileStartEntryIterator> & l,
                       const String<ChromosomeProfile::TChromosomeProfileStartEntryIterator> & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (position(l[i]) != position(r[i]))
            return false;
        if (!atEnd(l[i]) && !atEnd(r[i]))
        {
            if (!(*(l[i]) == *(r[i])))
            {
                return false;
            }
        }
        else if (atEnd(l[i]) ^ atEnd(r[i])) // Exactly one of them is atEnd
            return false;
    }
    return true;
}
inline bool operator==(const String<ChromosomeProfile::TChromosomeProfileEndEntryIterator> & l,
                       const String<ChromosomeProfile::TChromosomeProfileEndEntryIterator> & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (position(l[i]) != position(r[i]))
            return false;
        if (!atEnd(l[i]) && !atEnd(r[i]))
        {
            if (!(*(l[i]) == *(r[i])))
            {
                return false;
            }
        }
        else if (atEnd(l[i]) ^ atEnd(r[i])) // Exactly one of them is atEnd
            return false;
    }
    return true;
}
inline bool operator==(const ChromosomeProfile::TActiveSet & l, const ChromosomeProfile::TActiveSet & r)
{
    if (l.size() != r.size())
        return false;
    ChromosomeProfile::TActiveSet::const_iterator itL(l.begin());
    ChromosomeProfile::TActiveSet::const_iterator itEnd(l.end());
    while (itL != itEnd)
    {
        if (!(r.count(*itL)))
            return false;
        ++itL;
    }
    return true;
}
inline bool operator==(const String<ChromosomeProfile::TActiveSet>  & l,
                       const String<ChromosomeProfile::TActiveSet>  & r)
{
    if (length(l) != length(r))
        return false;
    for (unsigned i = 0; i < length(l); ++i)
    {
        if (!(l[i] == r[i]))
            return false;
    }
    return true;
}
inline bool operator==(const ChromosomeProfile & l, const ChromosomeProfile & r)
{
        if (l.currentPos != r.currentPos)
            return false;
        else if (l.globalMinPos != r.globalMinPos)
            return false;
        else if (l.avgNewReadsPerWindow != r.avgNewReadsPerWindow)
            return false;
        else if (!(l.startProfiles == r.startProfiles))
            return false;
        else if (!(l.endProfiles == r.endProfiles))
            return false;
        else if (!(l.activeReads == r.activeReads))
            return false;
        else
            return true;
}
// =======================================================================================
// Function performSwitches()
// =======================================================================================
// Perform the necessary switches for all readGroups for advancing to the next subset.
inline void performSwitches(ChromosomeProfile & profile,
                            const TReadGroupIndices & rg)
{
    for (unsigned r = 0; r < length(rg); ++r)
    {
        unsigned i = rg[r];
        while (!atEnd(profile.endProfiles[i]))
            profile.nextEndWindow(i);

        profile.startProfiles[i].switchBothSets(profile.currentRightBorder);
        profile.endProfiles[i].switchWriteSet();
        profile.endProfiles[i].correctConsecutiveSwitch(profile.activeLoad[i]);
        profile.endProfiles[i].switchReadSet();
    }
}
// =======================================================================================
// Function performPartialSwitches()
// =======================================================================================
// Perform the necessary switches for all readGroups for advancing to the next subset.
// But do NOT switch the writeSet of the end profiles.
inline void performPartialSwitches(ChromosomeProfile & profile,
                                  const TReadGroupIndices & rg)
{
    for (unsigned r = 0; r < length(rg); ++r)
    {
        unsigned i = rg[r];
        while (!atEnd(profile.endProfiles[i]))
            profile.nextEndWindow(i);

        profile.startProfiles[i].switchBothSets(profile.currentRightBorder);
        profile.endProfiles[i].switchReadSet();
    }
}
// =======================================================================================
// Function checkAllEmpty()
// =======================================================================================
// Return true if all start and end tables for all read groups are empty, false otherwise.
inline bool checkAllEmpty(const ChromosomeProfile & profile)
{
    bool allEmpty = true;
    for (unsigned i = 0; i < length(profile.startProfiles); ++i)
    {
        allEmpty &= empty(profile.startProfiles[i].sets[profile.startProfiles[i].writeSet].subset);
        allEmpty &= empty(profile.endProfiles[i].sets[profile.endProfiles[i].writeSet].subset);
        allEmpty &= profile.activeReads[i].empty();
        if (!allEmpty)
            return false;
    }
    return true;
}
#endif /* PROFILE_STRUCTURE_POPDEL_H_ */
