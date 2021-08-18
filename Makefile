CC      =gcc
CFLAGS  =-O3 -Wall 
LDFLAGS =-lm
SRCDIR  =src
OBJDIR  =obj
TARGET1 =example1.out
TRGSRC1 =example1.c
TARGET2 =example2.out
TRGSRC2 =example2.c

SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS)) ))
HEAD=$(wildcard $(SRCDIR)/*.h)

TRGOBJ1=$(OBJS) $(patsubst %.c,%.o,$(TRGSRC1))
TRGOBJ2=$(OBJS) $(patsubst %.c,%.o,$(TRGSRC2))


all : directories $(TARGET1) $(TARGET2)

directories:
	@mkdir -p $(OBJDIR)

$(TARGET1) : $(TRGOBJ1) 
	$(CC) $(LDFLAGS) -o $@ $^

$(TARGET2) : $(TRGOBJ2) 
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I $(SRCDIR) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) -I$(SRCDIR) -c $<

clean:
	@rm -rf $(TARGET1) $(TARGET2) $(OBJDIR) ./*.o

$(OBJS) : $(HEAD)
