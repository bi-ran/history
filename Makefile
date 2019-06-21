CXX = g++
CXXFLAGS += -O2 -Wall -Werror -Wextra -std=c++14
ROOTFLAGS := `root-config --cflags --libs`

BLDDIR = ./build
LIBDIR = ./lib
SRCDIR = ./src

LIBCONF = libhist.a

SRCS = $(wildcard $(SRCDIR)/*.C)
OBJS = $(patsubst $(SRCDIR)/%.C,$(BLDDIR)/%.o,$(SRCS))
DEPS = $(patsubst $(SRCDIR)/%.C,$(BLDDIR)/%.d,$(SRCS))

.PHONY: clean

$(LIBDIR)/$(LIBCONF): $(OBJS)
	@mkdir -p $(LIBDIR)
	ar rc $@ $^
	ranlib $@

$(BLDDIR)/%.o : $(SRCDIR)/%.C
	@mkdir -p $(BLDDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -MMD -MF $(BLDDIR)/$(*F).d $< -c -o $@

clean:
	@$(RM) $(LIBDIR)/$(LIBCONF) $(OBJS) $(DEPS)
	@rm -rf $(BLDDIR)/*

-include $(DEPS)
