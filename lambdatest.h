// *************************************************************************************************
// The MIT License (MIT)
// 
// Copyright (c) 2014 Christopher James Huff
//                    https://github.com/cjameshuff/
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// *************************************************************************************************

//! \file simpletest.h
//! \brief A simple unit testing framework for C++11+, making use of templates and lambdas.
//! Testing is organized as test groups, tests, and assertions.
//!
//! Test groups may be nested, and support setup and teardown handlers that are executed before and
//! after every test. These handlers may be redefined between the contained tests, or nested test
//! groups used to restrict them to a subset of tests.
//!
//! TODO: Better reporting reflecting group/test structure

#include <string>
#include <cstdio>
#include <chrono>
#include <functional>
#include <list>
#include <vector>

namespace lambdatest {

#define LAMBDATEST_DECL  namespace lambdatest {TestManager lambdatest_testManager;}

// *************************************************************************************************
// Display attributes
//      0 Reset all attributes
//      1 Bright
//      2 Dim
//      4 Underscore 
//      5 Blink
//      7 Reverse
//      8 Hidden
// Foreground colors:
//     30 Black
//     31 Red
//     32 Green
//     33 Yellow
//     34 Blue
//     35 Magenta
//     36 Cyan
//     37 White
// Background Colours
//     40 Black
//     41 Red
//     42 Green
//     43 Yellow
//     44 Blue
//     45 Magenta
//     46 Cyan
//     47 White (gray)

#define FAIL_TEST_ATTR  "\033[1;31;43m"
#define PASS_TEST_ATTR  "\033[1;32m"
#define FAIL_SUMM_ATTR  "\033[1;31m"
#define PASS_SUMM_ATTR  "\033[1;32m"
#define RESET_ATTR  "\033[0m"

struct TestContext {
    std::string groupName;
    bool stop;
    bool stopOnFail;
    std::function<void()> setup;
    std::function<void()> teardown;
};

struct LogEntry {
    std::string groupName;
    int testNumber;
    int assertionNumber;
    std::string description;
    bool passed;
};

struct TestManager {
    int testCounter = 0;
    int assertionCounter = 0;
    int assertionsPassed = 0;
    int assertionsFailed = 0;
    
    std::vector<TestContext> testContexts;
    std::list<LogEntry> testLog;
    
    auto SummarizeTests() -> void;
    
    TestManager() {}
    ~TestManager() {SummarizeTests();}
    
    auto SetupTest() -> void;
    auto TeardownTest() -> void;
    
    template<typename T>
    auto TestGroup(const std::string & groupName, T tests) -> void;
    
    template<typename Test>
    auto Test(const std::string & testName, Test test) -> void;
    
    auto ReportAssertionResult(const std::string & desc, bool passed) -> void;
};

extern TestManager lambdatest_testManager;


inline auto CurrentTestGroup() -> TestContext & {return lambdatest_testManager.testContexts.back();}

inline auto StopOnFail(bool sof) -> void {CurrentTestGroup().stopOnFail = sof;}

template<typename Handler>
inline auto Setup(Handler sh) -> void {CurrentTestGroup().setup = sh;}

template<typename Handler>
inline auto Teardown(Handler ch) -> void {CurrentTestGroup().teardown = ch;}


//! \brief Call setup handlers for test.
inline auto TestManager::SetupTest() -> void
{
    for(auto su = testContexts.begin(); su != testContexts.end(); ++su)
    {
        try {
            if(su->setup)
                su->setup();
        }
        catch(std::exception err) {
            std::clog << "Unexpected exception caught in test setup for test group \""
                      << su->groupName <<  "\"." << std::endl;
            std::clog << err.what() << std::endl;
            throw;
        }
        catch(...) {
            std::clog << "Unexpected exception caught in test setup for test group \""
                      << su->groupName <<  "\"." << std::endl;
            throw;
        }
    }
}

//! \brief Call teardown handlers for test.
inline auto TestManager::TeardownTest() -> void
{
    for(auto su = testContexts.rbegin(); su != testContexts.rend(); ++su)
    {
        try {
            if(su->teardown)
                su->teardown();
        }
        catch(std::exception err) {
            std::clog << "Unexpected exception caught in test teardown for test group \""
                      << su->groupName <<  "\"." << std::endl;
            std::clog << err.what() << std::endl;
            throw;
        }
        catch(...) {
            std::clog << "Unexpected exception caught in test teardown for test group \""
                      << su->groupName <<  "\"." << std::endl;
            throw;
        }
    }
}


//! \brief Perform a test, first running all setup tasks and running all cleanup tasks afterward.
//! If an exception is encountered during test setup or teardown, testing is aborted immediately.
template<typename TestFn>
auto TestManager::Test(const std::string & testName, TestFn test) -> void
{
    std::clog << "\033[1;34m""Test " << testCounter << " \"" << testName << "\":";
    std::clog << RESET_ATTR << std::endl;
    auto startT = std::chrono::steady_clock::now();
    
    assertionCounter = 0;
    SetupTest();
    
    test();
    
    TeardownTest();
    
    auto elapsed = std::chrono::duration<double>{std::chrono::steady_clock::now() - startT};
    std::clog << "\033[1;34m""Test " << testCounter << " \"" << testName
              << "\" completed in " << elapsed.count() << " s.";
    std::clog << RESET_ATTR << std::endl << std::endl;
    
    ++testCounter;
}

template<typename TestFn>
auto Test(const std::string & testName, TestFn test) -> void
{
    lambdatest_testManager.Test(testName, test);
}


//! \brief Define a test group.
template<typename T>
auto TestManager::TestGroup(const std::string & groupName, T tests) -> void
{
    auto extendedGroupName = groupName;
    if(!testContexts.empty())
    {
        extendedGroupName = CurrentTestGroup().groupName + ":" + groupName;
    }
    
    testContexts.push_back({extendedGroupName, false, false});
    
    std::clog << "********************************************************************************"
              << std::endl;
    std::clog << "Test group \"" << groupName << "\":" << std::endl;
    std::clog << "--------------------------------------------------------------------------------"
              << std::endl;
    
    auto startT = std::chrono::steady_clock::now();
    
    try {
        tests();
    }
    catch(std::exception err) {
        std::clog << "Unexpected exception caught in TestGroup() " << groupName << ":" << std::endl;
        std::clog << err.what() << std::endl;
        throw;
    }
    catch(...) {
        std::clog << "Unexpected exception caught in TestGroup()" << groupName << ":" << std::endl;
        throw;
    }
    
    auto elapsed = std::chrono::duration<double>{std::chrono::steady_clock::now() - startT};
    
    std::clog << "--------------------------------------------------------------------------------"
              << std::endl;
    std::clog << "Test group \"" << groupName
              << "\" completed in " << elapsed.count() << " s."
              << std::endl;
    
    std::clog << "********************************************************************************"
              << std::endl << std::endl;
    
    testContexts.pop_back();
}

template<typename T>
auto TestGroup(const std::string & groupName, T tests) -> void
{
    lambdatest_testManager.TestGroup(groupName, tests);
}


//! \brief Report result of an assertion, with string giving a brief description.
inline auto TestManager::ReportAssertionResult(const std::string & desc, bool passed) -> void
{
    testLog.push_back({
        CurrentTestGroup().groupName,
        testCounter,
        assertionCounter,
        desc,
        passed
    });
    
    if(passed)
    {
        std::clog << PASS_TEST_ATTR << "PASS: ";
        ++assertionsPassed;
    }
    else
    {
        std::clog << FAIL_TEST_ATTR << "FAIL: ";
        ++assertionsFailed;
        CurrentTestGroup().stop = CurrentTestGroup().stopOnFail;
    }
    std::clog << CurrentTestGroup().groupName
              << "[" << testCounter << ":" << assertionCounter << "] " << desc;
    std::clog << RESET_ATTR << std::endl;
    
    ++assertionCounter;
}


//! \brief Report summary of test results, with pass/fail counts.
inline auto TestManager::SummarizeTests() -> void
{
    std::clog << "\n----------------------------------------" << std::endl;
    std::clog << "Summary:" << std::endl;
    std::clog << "----------------------------------------" << std::endl;
    for(auto & le: testLog)
    {
        if(le.passed)
        {
            std::clog << PASS_TEST_ATTR;
            std::clog << "PASS: ";
        }
        else
        {
            std::clog << FAIL_TEST_ATTR;
            std::clog << "FAIL: ";
        }
        std::clog << le.groupName
                  << "[" << le.testNumber << ":" << le.assertionNumber << "] "
                  << le.description;
        std::clog << RESET_ATTR << std::endl;
    }
    std::clog << "----------------------------------------" << std::endl;
    std::clog << ((assertionsFailed)? FAIL_SUMM_ATTR : PASS_SUMM_ATTR);
    std::clog << assertionsPassed << " assertions passed, "
              << assertionsFailed << " assertions failed" << std::endl;
    std::clog << RESET_ATTR;
}


// *************************************************************************************************

//! \brief Assert that some operation should produce a true result.
template<typename Test>
inline auto Should(const std::string & desc, const Test & test) -> void {
    if(CurrentTestGroup().stop)
        return;
    
    bool passed = false;
    try {passed = test();}
    catch(std::exception err) {
        std::clog << "Unexpected exception caught" << std::endl;
        std::clog << err.what() << std::endl;
        passed = false;
    }
    catch(...) {
        std::clog << "Unexpected exception caught" << std::endl;
        passed = false;
    }
    lambdatest_testManager.ReportAssertionResult(desc, passed);
}

inline auto Should(const std::string & desc, bool test) -> void {
    Should(desc, [=]{return test;});
}

// *************************************************************************************************

template<typename T>
inline auto default_formatter(T value) -> void
{
    std::clog << value;
}

// *************************************************************************************************

//! \brief Variation of ShouldEq() assertion for an array of multiple values with custom formatter.
template<typename T, typename Formatter>
inline auto MShouldEq(std::string desc, T * value, T * expected, size_t count,
                      Formatter fmt = default_formatter) -> void
{
    Should(desc, [&]{
        bool pass = true;
        for(size_t j = 0; j < count && pass; ++j)
            pass = value[j] == expected[j];
        
        if(!pass)
        {
            for(size_t j = 0; j < count; ++j)
            {
                bool mismatch = value[j] != expected[j];
                std::clog << (mismatch? FAIL_TEST_ATTR: PASS_TEST_ATTR);
                std::clog << j << "  ";
                fmt(value[j]);
                std::clog << (mismatch? " != ": " == ");
                fmt(expected[j]);
                std::clog << RESET_ATTR;
                std::clog << std::endl;
            }
        }
        
        return pass;
    });
}

//! \brief Variation of ShouldEq() assertion for an array of multiple values.
template<typename T>
inline auto MShouldEq(std::string desc, T * value, T * expected, size_t count) -> void
{
    MShouldEq(desc, value, expected, count, [](T x){std::clog << x;});
}


// *************************************************************************************************
//! \brief Should assertion for binary operations.
template<typename T, typename T2, typename Pred, typename Formatter>
inline auto ShouldBin(std::string desc, T value, T2 expected,
                      Pred pred, std::string predStr, Formatter fmt) -> void
{
    Should(desc, [&]{
        bool pass = pred(value, expected);
        
        if(!pass)
        {
            std::clog << FAIL_TEST_ATTR << "Failed assertion: ";
            fmt(value);
            std::clog << predStr;
            fmt(expected);
            std::clog << RESET_ATTR << std::endl;
        }
        
        return pass;
    });
}


// *************************************************************************************************
//! \brief Should assertion comparing two values for equality.
template<typename T, typename T2, typename Formatter>
inline auto ShouldEq(std::string desc, T value, T2 expected, Formatter fmt) -> void
{
    ShouldBin(desc, value, expected, std::equal_to<T>(), " == ", fmt);
}

template<typename T, typename T2>
inline auto ShouldEq(std::string desc, T value, T2 expected) -> void
{
    ShouldBin(desc, value, expected, std::equal_to<T>(), " == ", default_formatter<T>);
}


// *************************************************************************************************
//! \brief Should assertion comparing two values for inequality.
template<typename T, typename T2, typename Formatter>
inline auto ShouldNEq(std::string desc, T value, T2 expected, Formatter fmt) -> void
{
    ShouldBin(desc, value, expected, std::not_equal_to<T>(), " != ", fmt);
}

template<typename T, typename T2>
inline auto ShouldNEq(std::string desc, T value, T2 expected) -> void
{
    ShouldBin(desc, value, expected, std::not_equal_to<T>(), " != ", default_formatter<T>);
}


// *************************************************************************************************

//! \brief Assert that code should not throw any exception.
template<typename Test>
auto ShouldNotThrow(std::string desc, const Test & test) -> void {
    if(CurrentTestGroup().stop)
        return;
    
    bool passed = true;
    try {test();}
    catch(...) {
        std::clog << "Unexpected exception caught" << std::endl;
        passed = false;
    }
    lambdatest_testManager.ReportAssertionResult(desc, passed);
}


//! \brief Assert that code should throw any exception.
template<typename Test>
auto ShouldThrowAny(const std::string & desc, const Test & test) -> void {
    if(CurrentTestGroup().stop)
        return;
    
    bool passed = false;
    try {test();}
    catch(...) {
        std::clog << "Caught expected exception" << std::endl;
        passed = true;
    }
    lambdatest_testManager.ReportAssertionResult(desc, passed);
}


//! \brief Assert that code should not throw a specific exception.
template<typename Excep, typename Test>
auto ShouldThrow(const std::string & desc, const Test & test) -> void {
    if(CurrentTestGroup().stop)
        return;
    
    bool passed = false;
    try {test();}
    catch(Excep & excep) {
        std::clog << "Caught expected exception" << std::endl;
        passed = true;
    }
    catch(...) {
        std::clog << "Unexpected exception caught" << std::endl;
        passed = false;
    }
    lambdatest_testManager.ReportAssertionResult(desc, passed);
}

// *************************************************************************************************
} // namespace lambdatest
