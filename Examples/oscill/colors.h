namespace Gascoigne
{
namespace style
{
/*! \brief Some escape sequences for output styling
 * Normally works on linux machines, won't work on windows machines. Terminal needs to understand
 * ANSI escape sequences.
 */
constexpr auto r = "\033[31m";  //!< red
constexpr auto g = "\033[32m";  //!< green
constexpr auto b = "\033[34m";  //!< blue
constexpr auto c = "\033[36m";  //!< cyan
constexpr auto m = "\033[35m";  //!< magenta
constexpr auto y = "\033[33m";  //!< yellow
constexpr auto k = "\033[30m";  //!< black
constexpr auto w = "\033[37m";  //!< white

constexpr auto R = "\033[91m";  //!< hi-red
constexpr auto G = "\033[92m";  //!< hi-green
constexpr auto B = "\033[94m";  //!< hi-blue
constexpr auto C = "\033[96m";  //!< hi-cyan
constexpr auto M = "\033[95m";  //!< hi-magenta
constexpr auto Y = "\033[99m";  //!< hi-yellow
constexpr auto K = "\033[90m";  //!< hi-black
constexpr auto W = "\033[97m";  //!< hi-white

constexpr auto bb = "\033[1m";   //!< bold
constexpr auto ul = "\033[4m";   //!< underline
constexpr auto ft = "\033[2m";   //!< faint
constexpr auto i  = "\033[7m";   //!< inverse
constexpr auto p  = "\033[27m";  //!< positive
constexpr auto n  = "\033[0m";   //!< normal
}  // namespace style
}  // namespace Gascoigne
