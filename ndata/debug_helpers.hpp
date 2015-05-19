#include <iostream>
#include <sstream>
#include <istream>

#include <tuple>

using namespace std;

#define DECLARE_TESTRESULT(acc_success_status, acc_message) \
    bool acc_success_status=true; \
    string acc_message (""); \
    bool success_status_TMP =true; \
    acc_success_status = success_status_TMP; /*disable warning when not using TEST macro*/ \
    string messages_TMP (""); \
    acc_message = messages_TMP; \
    \

#define TEST(test_func, acc_success_status, acc_message) \
    acc_message.append(#test_func "\t"); \
    tie(success_status_TMP, messages_TMP) = test_func; \
    acc_success_status = acc_success_status and success_status_TMP; \
    if (success_status_TMP==true) { \
        acc_message.append("success\n"); \
    } else { \
        acc_message.append("FAILURE"); \
    } \
    /* add message from test to acc_message*/ \
    acc_message.append(messages_TMP); \
    acc_message.append("\n"); \
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
 * Use like that
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


typedef pair<bool, string> TestResult;
