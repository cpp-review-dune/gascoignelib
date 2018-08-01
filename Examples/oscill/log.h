namespace Gascoigne
{
namespace logwrite
{
//! \brief Some helper functions for simplifying log writing
namespace implementation
{
//! \brief Internal functions for logwriting
template <typename T>
std::string toString(T value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

std::string merge(std::initializer_list<std::string> strList)
{
    std::string ret = "";
    for (std::string s : strList)
    {
        ret += s;
    }
    return ret;
}
}  // namespace implementation

using namespace implementation;

template <typename... slist>
std::string stringer(const slist&... tbstr)
{
    return merge({toString(tbstr)...});
}

template <typename... slist>
std::string info(const std::string& info_name, const slist&... messages)
{
    return merge({toString(sty::bb), "[", info_name, "] ", toString(sty::n),
                  toString(messages)..., "\n", toString(sty::n), toString(sty::k)});
}

}  // namespace logwrite
}  // namespace Gascoigne
