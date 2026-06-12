/**
 * PRadTrigger body file. That allows for decoding of the PRadTrigger vector.
 *
 * @author - Rafayel Paramuzyan <rafo@jlab.org>
 * @see - https://code.jlab.org/hallb/prad_x17/pradx17_trigger/-/blob/main/src/PRadTrigger.cc?ref_type=heads
 */

#include "PRadTrigger.h"

#include <iostream>
#include <TTree.h>

// nglobal_trigger_words is stored as UChar_t (/b) in the PRad2 ROOT files
// → maximum 255 words per event
static constexpr int kMaxSSPWords = 256;

// -----------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------

PRadTrigger::PRadTrigger()
    : m_eventNumber(-1)
    , m_triggerTimeNs(-1LL)
{
    // m_activeTimes and m_bitTimes[] are default-constructed (empty)
}

// -----------------------------------------------------------------------
// Reset
// -----------------------------------------------------------------------

void PRadTrigger::Clear()
{
    m_eventNumber   = -1;
    m_triggerTimeNs = -1LL;
    m_activeTimes.clear();
    for (int b = 0; b < N_TRIG_BITS; ++b)
        m_bitTimes[b].clear();
}

// -----------------------------------------------------------------------
// Parsing
// -----------------------------------------------------------------------

void PRadTrigger::Parse(const uint32_t* words, int nwords)
{
    Clear();
    if (!words || nwords <= 0) return;

    for (int i = 0; i < nwords; /* advanced inside */) {

        const uint32_t w0   = words[i];
        const uint32_t type = w0 & TYPE_MASK;

        // --- EVTHDR -----------------------------------------------------
        if (type == EVTHDR_VAL) {
            m_eventNumber = static_cast<int>(w0 & EVENTNUM_MASK);
            ++i;
            continue;
        }

        // --- TRGTIME ----------------------------------------------------
        if (type == TRGTIME_VAL) {
            if (i + 1 >= nwords) {
                std::cerr << "[PRadTrigger] Truncated TRGTIME at word " << i << "\n";
                break;
            }
            const uint32_t w1    = words[i + 1];
            const uint64_t timeL = static_cast<uint64_t>(w0 & TRGTIME_LO_MASK);
            const uint64_t timeH = static_cast<uint64_t>(w1 & TRGTIME_HI_MASK);
            // 48-bit tick counter; each tick = 4 ns
            const uint64_t ticks = (timeH << 24) | timeL;
            m_triggerTimeNs = static_cast<long long>(ticks) * 4LL;
            i += 2;
            continue;
        }

        // --- TRIGGER ----------------------------------------------------
        if (type == TRIGGER_VAL) {
            const uint32_t trigTime = (w0 & TRIG_TIME_MASK) >> TRIG_TIME_SHIFT;
            const uint32_t pattern  =  w0 & TRIG_PAT_MASK;

            if (pattern != 0) {
                m_activeTimes.insert(trigTime);
                for (int bit = 0; bit < N_TRIG_BITS; ++bit) {
                    if (pattern & (1u << bit))
                        m_bitTimes[bit].push_back(static_cast<int>(trigTime));
                }
            }
            ++i; // single word only
            continue;
        }

        // --- anything else: skip ----------------------------------------
        ++i;
    }
}

void PRadTrigger::Parse(const std::vector<unsigned int>& words)
{
    Parse(reinterpret_cast<const uint32_t*>(words.data()),
          static_cast<int>(words.size()));
}

void PRadTrigger::Parse(const std::vector<int>& words)
{
    // Same bit layout as the unsigned variant; reinterpret.
    Parse(reinterpret_cast<const uint32_t*>(words.data()),
          static_cast<int>(words.size()));
}

// -----------------------------------------------------------------------
// Getters — EVTHDR
// -----------------------------------------------------------------------

int PRadTrigger::GetEventNumber() const
{
    return m_eventNumber;
}

// -----------------------------------------------------------------------
// Getters — TRGTIME
// -----------------------------------------------------------------------

long long PRadTrigger::GetEventTriggerTime() const
{
    return m_triggerTimeNs;
}

// -----------------------------------------------------------------------
// Getters — TRIGGER
// -----------------------------------------------------------------------

int PRadTrigger::GetNTrigTimes() const
{
    return static_cast<int>(m_activeTimes.size());
}

const std::set<uint32_t>& PRadTrigger::GetActiveTrigTimes() const
{
    return m_activeTimes;
}

const std::vector<int>* PRadTrigger::GetTriggersBit(int ibit) const
{
    if (ibit < 0 || ibit >= N_TRIG_BITS) return nullptr;
    return &m_bitTimes[ibit];
}

int PRadTrigger::GetNFirings(int ibit) const
{
    const auto* v = GetTriggersBit(ibit);
    return v ? static_cast<int>(v->size()) : 0;
}

bool PRadTrigger::IsBitActive(int ibit) const
{
    return GetNFirings(ibit) > 0;
}

// -----------------------------------------------------------------------
// Debug
// -----------------------------------------------------------------------

void PRadTrigger::Print() const
{
    std::cout << "=== PRadTrigger ===\n"
              << "  EventNumber       : " << m_eventNumber   << "\n"
              << "  EventTrigTime(ns) : " << m_triggerTimeNs << "\n"
              << "  NTrigTimes        : " << GetNTrigTimes() << "\n";
    for (int b = 0; b < N_TRIG_BITS; ++b) {
        if (!m_bitTimes[b].empty()) {
            std::cout << "  Bit " << b << " -> times: ";
            for (int t : m_bitTimes[b]) std::cout << t << " ";
            std::cout << "\n";
        }
    }
}
