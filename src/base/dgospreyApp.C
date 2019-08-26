#include "dgospreyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<dgospreyApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

dgospreyApp::dgospreyApp(InputParameters parameters) : MooseApp(parameters)
{
  dgospreyApp::registerAll(_factory, _action_factory, _syntax);
}

dgospreyApp::~dgospreyApp() {}

void
dgospreyApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"dgospreyApp"});
  Registry::registerActionsTo(af, {"dgospreyApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
dgospreyApp::registerApps()
{
  registerApp(dgospreyApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
dgospreyApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  dgospreyApp::registerAll(f, af, s);
}
extern "C" void
dgospreyApp__registerApps()
{
  dgospreyApp::registerApps();
}
