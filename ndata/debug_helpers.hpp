#include <iostream>
#include <sstream>
#include <istream>

#include <tuple>

using namespace std;

//comment/uncomment to control printing behavior of the test macro
//#define PRINT_TESTRESULT_EVEN_ON_SUCCESS

#ifdef PRINT_TESTRESULT_EVEN_ON_SUCCESS
#define MAYBE_APPEND_MESSAGE_ON_FAILURE(target, toappend)
#define MAYBE_APPEND_MESSAGE_ALWAYS(target, toappend) target.append(toappend)
#else
#define MAYBE_APPEND_MESSAGE_ON_FAILURE(target, toappend) target.append(toappend)
#define MAYBE_APPEND_MESSAGE_ALWAYS(target, toappend)
#endif

#define DECLARE_TEST(acc_success_status, acc_message) \
    bool acc_success_status=true; \
    string acc_message (""); \
    bool success_status_TMP =true; \
    acc_success_status = success_status_TMP; /*disable unused variable warning when not using RUN_TEST macro*/ \
    string messages_TMP (""); \
    acc_message = messages_TMP; \
    \

#define RUN_TEST(test_func, acc_success_status, acc_message) \
    acc_message.append(#test_func "\t"); \
    tie(success_status_TMP, messages_TMP) = test_func; \
    acc_success_status = acc_success_status and success_status_TMP; \
    if (success_status_TMP==true) { \
        acc_message.append("success\n"); \
    } else { \
        acc_message.append("FAILURE\n"); \
        MAYBE_APPEND_MESSAGE_ON_FAILURE(acc_message, messages_TMP); \
    } \
    MAYBE_APPEND_MESSAGE_ALWAYS(acc_message, messages_TMP); \
    acc_message.append("\n"); \
    \

/**
 * nicely indents the result to reflect test function structure
 */
#define RETURN_TESTRESULT(acc_success_status, acc_message) \
    string ret_message (""); \
    istringstream acc_message_istream (acc_message); \
    for(string line; getline(acc_message_istream, line); ) { \
        ret_message.append(MakeString() << "\t" << line << "\n"); \
    } \
    return make_pair(acc_success_status, ret_message); \
    \

/**
 *
 * Use like that :
 *
 * #include <string>
 * #include <sstream>
 * #include <iostream>
 *
 * MakeString() << val << "stuff " << endl
 */
class MakeString {
    public:
        std::stringstream stream;
        operator std::string() const { return stream.str(); }

        template<class T>
        MakeString& operator<<(T const& VAR) { stream << VAR; return *this; }
};


typedef pair<bool, string> test_result;
