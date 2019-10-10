include Make.inc

all : $(OBJS)
	$(CXX) $(CXXFLAGS) -l$(LIBS) -o $(GOAL) $(OBJS)





%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I $(INCLUDE_PATH)

.PHONY: clean
clean:
	$(RM) $(OBJS)
	$(RM) $(GOAL)
