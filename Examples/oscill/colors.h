namespace Gascoigne::style
{
/*! \brief Some escape sequences for output styling
 * Normally works on linux machines, won't work on windows machines. Terminal needs to understand
 * ANSI escape sequences.
 */
const char* const r = "\033[31m";  //!< red
const char* const g = "\033[32m";  //!< green
const char* const b = "\033[34m";  //!< blue
const char* const c = "\033[36m";  //!< cyan
const char* const m = "\033[35m";  //!< magenta
const char* const y = "\033[33m";  //!< yellow
const char* const k = "\033[30m";  //!< black
const char* const w = "\033[37m";  //!< white

const char* const R = "\033[91m";  //!< hi-red
const char* const G = "\033[92m";  //!< hi-green
const char* const B = "\033[94m";  //!< hi-blue
const char* const C = "\033[96m";  //!< hi-cyan
const char* const M = "\033[95m";  //!< hi-magenta
const char* const Y = "\033[99m";  //!< hi-yellow
const char* const K = "\033[90m";  //!< hi-black
const char* const W = "\033[97m";  //!< hi-white

const char* const bb = "\033[1m";   //!< bold
const char* const ul = "\033[4m";   //!< underline
const char* const ft = "\033[2m";   //!< faint
const char* const i  = "\033[7m";   //!< inverse
const char* const p  = "\033[27m";  //!< positive
const char* const n  = "\033[0m";   //!< normal
}  // namespace Gascoigne::style
