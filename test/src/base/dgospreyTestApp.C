//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "dgospreyTestApp.h"
#include "dgospreyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<dgospreyTestApp>()
{
  InputParameters params = validParams<dgospreyApp>();
  return params;
}

dgospreyTestApp::dgospreyTestApp(InputParameters parameters) : MooseApp(parameters)
{
  dgospreyTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

dgospreyTestApp::~dgospreyTestApp() {}

void
dgospreyTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  dgospreyApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"dgospreyTestApp"});
    Registry::registerActionsTo(af, {"dgospreyTestApp"});
  }
}

void
dgospreyTestApp::registerApps()
{
  registerApp(dgospreyApp);
  registerApp(dgospreyTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
dgospreyTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  dgospreyTestApp::registerAll(f, af, s);
}
extern "C" void
dgospreyTestApp__registerApps()
{
  dgospreyTestApp::registerApps();
}
