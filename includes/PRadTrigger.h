/**
 * PRadTrigger header file. That allows for decoding of the PRadTrigger vector.
 *
 * @author - Rafayel Paramuzyan <rafo@jlab.org>
 * @see - https://code.jlab.org/hallb/prad_x17/pradx17_trigger/-/blob/main/include/PRadTrigger.h?ref_type=heads
 */

#ifndef PRADTRIGGER_H
#define PRADTRIGGER_H

#include <vector>
#include <set>
#include <cstdint>

class TTree;

/**
 * PRadTrigger
 * Parses SSP Hardware Data words from the global_trigger_raw array.
 *
 * Handled word types:
 * -----------------------------------------------------------------------
 * EVTHDR  0(31:27) = 0x10+0x02 = 0x12  ->  word & TYPE_MASK = 0x90000000
 *         bits[26:0]  = EVENTNUM
 *
 * TRGTIME 0(31:27) = 0x10+0x03 = 0x13  ->  word & TYPE_MASK = 0x98000000
 *         word0 bits[23:0] = TIMEL
 *         word1 bits[23:0] = TIMEH
 *         full time = (uint64_t(TIMEH) << 24) | TIMEL  [1 tick = 4 ns]
 *
 * TRIGGER 0(31:27) = 0x10+0x0D = 0x1D  ->  word & TYPE_MASK = 0xE8000000
 *         bits[26:16] = TRIG_TIME     (local SSP trigger time, 11-bit)
 *         bits[15:0]  = TRIG_PATTERN  (single word; multiple bits may be set
 *                                      simultaneously for a given TRIG_TIME)
 * -----------------------------------------------------------------------
 *
 * PRad2 ROOT file branch notes:
 *   Old format:
 *     "nglobal_trigger_words" -> UChar_t  (/b)  max 255 words per event
 *     "global_trigger_raw"    -> UInt_t[] (/i)
 *   New format:
 *     "ssp_raw"               -> std::vector<unsigned int>
 *                                (legacy variant: std::vector<int> — same bit layout)
 */
class PRadTrigger {
public:

    // -----------------------------------------------------------------------
    // Word-type masks & values  (bits 31:27)
    // -----------------------------------------------------------------------
    static constexpr uint32_t TYPE_MASK   = 0xF8000000u;  ///< bits[31:27] selector

    static constexpr uint32_t EVTHDR_VAL  = 0x90000000u;  ///< 0x12 << 27
    static constexpr uint32_t TRGTIME_VAL = 0x98000000u;  ///< 0x13 << 27
    static constexpr uint32_t TRIGGER_VAL = 0xE8000000u;  ///< 0x1D << 27

    // EVTHDR field
    static constexpr uint32_t EVENTNUM_MASK   = 0x07FFFFFFu; ///< bits[26:0]

    // TRGTIME fields
    static constexpr uint32_t TRGTIME_LO_MASK = 0x00FFFFFFu; ///< bits[23:0] of word0
    static constexpr uint32_t TRGTIME_HI_MASK = 0x00FFFFFFu; ///< bits[23:0] of word1

    // TRIGGER fields
    static constexpr uint32_t TRIG_TIME_MASK  = 0x07FF0000u; ///< bits[26:16]
    static constexpr int      TRIG_TIME_SHIFT = 16;
    static constexpr uint32_t TRIG_PAT_MASK   = 0x0000FFFFu; ///< bits[15:0]
    static constexpr int      N_TRIG_BITS     = 16;           ///< pattern is 16-bit

    // -----------------------------------------------------------------------
    // Construction / destruction
    // -----------------------------------------------------------------------
    PRadTrigger();
    ~PRadTrigger() = default;

    // -----------------------------------------------------------------------
    // Reset
    // -----------------------------------------------------------------------
    /// Clear all parsed state. Called automatically at the start of Parse().
    void Clear();

    // -----------------------------------------------------------------------
    // Parsing
    // -----------------------------------------------------------------------

    /**
     * Parse one event's worth of raw SSP words (old array format).
     *
     * @param words   Pointer to the global_trigger_raw array
     * @param nwords  Value of nglobal_trigger_words (cast from UChar_t)
     */
    void Parse(const uint32_t* words, int nwords);

    /**
     * Parse one event's worth of raw SSP words (new vector format).
     *
     * @param words  Reference to the ssp_raw vector<unsigned int> branch
     */
    void Parse(const std::vector<unsigned int>& words);

    /// Overload for vector<int>-typed ssp_raw branches (same bit layout).
    void Parse(const std::vector<int>& words);

    // -----------------------------------------------------------------------
    // Getters — EVTHDR
    // -----------------------------------------------------------------------

    /**
     * Event number decoded from the EVTHDR word (bits[26:0]).
     * @return Event number, or -1 if no EVTHDR word was present.
     */
    int GetEventNumber() const;

    // -----------------------------------------------------------------------
    // Getters — TRGTIME
    // -----------------------------------------------------------------------

    /**
     * Global event trigger time in nanoseconds, reconstructed from the
     * 48-bit TRGTIME tick counter (TIMEH:TIMEL, 1 tick = 4 ns).
     * @return Time in ns, or -1 if no TRGTIME word pair was present.
     */
    long long GetEventTriggerTime() const;

    // -----------------------------------------------------------------------
    // Getters — TRIGGER
    // -----------------------------------------------------------------------

    /**
     * Number of distinct TRIG_TIME values that had at least one pattern
     * bit set.
     */
    int GetNTrigTimes() const;

    /**
     * Read-only access to the sorted set of all active TRIG_TIME values
     * (those for which the TRIG_PATTERN was non-zero).
     */
    const std::set<uint32_t>& GetActiveTrigTimes() const;

    /**
     * Returns a pointer to the internal vector of TRIG_TIME values at
     * which trigger bit @p ibit (0..15) was set.
     *
     * A single TRIGGER word can have multiple bits set simultaneously in
     * its TRIG_PATTERN field.  There is at most one TRIGGER word per
     * TRIG_TIME value, so the same TRIG_TIME can appear in the vectors
     * of several different bits, but never more than once in any single
     * bit's vector.
     *
     * @param ibit  Trigger bit index in [0, N_TRIG_BITS).
     * @return      Pointer to the time vector, or nullptr if @p ibit is
     *              out of range.
     */
    const std::vector<int>* GetTriggersBit(int ibit) const;

    /**
     * Number of trigger times at which bit @p ibit fired.
     * Returns 0 for an out-of-range @p ibit.
     */
    int GetNFirings(int ibit) const;

    /**
     * Returns true if trigger bit @p ibit fired at least once in the
     * last parsed event.
     */
    bool IsBitActive(int ibit) const;

    // -----------------------------------------------------------------------
    // Debug
    // -----------------------------------------------------------------------
    /// Print a compact summary of all parsed data to stdout.
    void Print() const;

private:
    // --- EVTHDR ---
    int       m_eventNumber;

    // --- TRGTIME (nanoseconds) ---
    long long m_triggerTimeNs;

    // --- TRIGGER ---
    std::set<uint32_t> m_activeTimes;           ///< sorted set of TRIG_TIME values with at least one bit set
    std::vector<int>   m_bitTimes[N_TRIG_BITS]; ///< per-bit list of TRIG_TIME values; one TRIG_TIME may appear across multiple bit vectors
};

#endif // PRADTRIGGER_H