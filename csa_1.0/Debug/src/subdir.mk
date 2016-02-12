################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BitRank.cpp \
../src/CHgtArray.cpp \
../src/CRMQ.cpp \
../src/CSA.cpp \
../src/Hash.cpp \
../src/LcpToParentheses.cpp \
../src/Parentheses.cpp \
../src/ReplacePattern.cpp \
../src/SSTree.cpp \
../src/SubblockRMQ.cpp \
../src/TestLCSS.cpp \
../src/Tools.cpp \
../src/bittree.cpp \
../src/rbtree.cpp \
../src/wtreebwt.cpp 

OBJS += \
./src/BitRank.o \
./src/CHgtArray.o \
./src/CRMQ.o \
./src/CSA.o \
./src/Hash.o \
./src/LcpToParentheses.o \
./src/Parentheses.o \
./src/ReplacePattern.o \
./src/SSTree.o \
./src/SubblockRMQ.o \
./src/TestLCSS.o \
./src/Tools.o \
./src/bittree.o \
./src/rbtree.o \
./src/wtreebwt.o 

CPP_DEPS += \
./src/BitRank.d \
./src/CHgtArray.d \
./src/CRMQ.d \
./src/CSA.d \
./src/Hash.d \
./src/LcpToParentheses.d \
./src/Parentheses.d \
./src/ReplacePattern.d \
./src/SSTree.d \
./src/SubblockRMQ.d \
./src/TestLCSS.d \
./src/Tools.d \
./src/bittree.d \
./src/rbtree.d \
./src/wtreebwt.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


