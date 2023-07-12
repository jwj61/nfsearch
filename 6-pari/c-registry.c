#include <strings.h>
#include "c-target.h"

registry_data registry[15];
static int reg_cnt=0;

void register_new_func(char *nm, int (*inf)(), long (*nf)()) {
  strcpy(registry[reg_cnt].name, nm);
  registry[reg_cnt].next_func = nf;
  registry[reg_cnt].init_func = inf;
  reg_cnt++;
}

void register_known_funcs() {
  register_new_func("t6", t6init, t6next);
  register_new_func("t51", t51init, t51next);
  register_new_func("t42", t42init, t42next);
  register_new_func("t411", t411init, t411next);
  register_new_func("t33", t33init, t33next);
  register_new_func("t321", t321init, t321next);
  register_new_func("t3111", t3111init, t3111next);
  register_new_func("t222", t222init, t222next);
  register_new_func("t2211", t2211init, t2211next);
  register_new_func("t21111", t21111init, t21111next);
}

int lookup_registry(char *nm) {
  int j=0;

  while(j<reg_cnt && !strcmp(nm, registry[j].name))
    j++;

  return (j<reg_cnt?j:-1);
}

  
  
