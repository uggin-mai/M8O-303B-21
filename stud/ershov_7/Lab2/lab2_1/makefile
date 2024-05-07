FLAGS_G++ = -g -pedantic -std=c++17 -Wall -Wextra -Werror
PROGRAN_NAME = solution


compile: main.cpp
	g++ $(FLAGS_G++) main.cpp -o $(PROGRAN_NAME)

run: compile
	$(PROGRAN_NAME) < input > output

clean:
	del $(PROGRAN_NAME).exe