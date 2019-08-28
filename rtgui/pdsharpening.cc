/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include "eventmapper.h"
#include "pdsharpening.h"
#include "options.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

PdSharpening::PdSharpening() : FoldableToolPanel(this, "pdsharpening", M("TP_PDSHARPENING_LABEL"), false, true)
{

    auto m = ProcEventMapper::getInstance();
    EvPdShrContrast = m->newEvent(ALLNORAW, "HISTORY_MSG_PDSHARPEN_CONTRAST");
    EvPdSharpenGamma = m->newEvent(ALLNORAW, "HISTORY_MSG_PDSHARPEN_GAMMA");
    EvPdShrDRadius = m->newEvent(ALLNORAW, "HISTORY_MSG_PDSHARPEN_RADIUS");
    EvPdShrDIterations = m->newEvent(ALLNORAW, "HISTORY_MSG_PDSHARPEN_ITERATIONS");
    EvPdShrAutoContrast = m->newEvent(ALLNORAW, "HISTORY_MSG_PDSHARPEN_AUTO_CONTRAST");

    Gtk::HBox* hb = Gtk::manage(new Gtk::HBox());
    hb->show();
    contrast = Gtk::manage(new Adjuster(M("TP_SHARPENING_CONTRAST"), 0, 200, 1, 10));
    contrast->setAdjusterListener(this);
    contrast->addAutoButton(M("TP_RAW_DUALDEMOSAICAUTOCONTRAST_TOOLTIP"));
    contrast->setAutoValue(true);

    pack_start(*contrast);
    contrast->show();

    pack_start(*hb);

    Gtk::VBox* rld = Gtk::manage(new Gtk::VBox());
    gamma = Gtk::manage(new Adjuster(M("TP_SHARPENING_GAMMA"), 0.5, 6.0, 0.05, 1.35));
    dradius = Gtk::manage(new Adjuster(M("TP_SHARPENING_EDRADIUS"), 0.4, 1.0, 0.01, 0.75));
    diter = Gtk::manage(new Adjuster(M("TP_SHARPENING_RLD_ITERATIONS"), 5, 100, 1, 30));
    rld->pack_start(*gamma);
    rld->pack_start(*dradius);
    rld->pack_start(*diter);
    gamma->show();
    dradius->show();
    diter->show();
    rld->show();
    pack_start(*rld);

    dradius->setAdjusterListener(this);
    gamma->setAdjusterListener(this);
    diter->setAdjusterListener(this);

    contrast->delay = std::max(contrast->delay, options.adjusterMaxDelay);
    dradius->delay = std::max(dradius->delay, options.adjusterMaxDelay);
    gamma->delay = std::max(gamma->delay, options.adjusterMaxDelay);
    diter->delay = std::max(diter->delay, options.adjusterMaxDelay);
}

PdSharpening::~PdSharpening()
{
    idle_register.destroy();
}


void PdSharpening::read(const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener();

    if (pedited) {
        contrast->setEditedState(pedited->pdsharpening.contrast ? Edited : UnEdited);
        contrast->setAutoInconsistent(multiImage && !pedited->pdsharpening.autoContrast);
        gamma->setEditedState(pedited->pdsharpening.gamma ? Edited : UnEdited);
        dradius->setEditedState(pedited->pdsharpening.deconvradius ? Edited : UnEdited);
        diter->setEditedState(pedited->pdsharpening.deconviter ? Edited : UnEdited);

        set_inconsistent(multiImage && !pedited->pdsharpening.enabled);
    }

    setEnabled(pp->pdsharpening.enabled);

    contrast->setValue(pp->pdsharpening.contrast);
    contrast->setAutoValue(pp->pdsharpening.autoContrast);
    gamma->setValue(pp->pdsharpening.gamma);
    dradius->setValue(pp->pdsharpening.deconvradius);
    diter->setValue(pp->pdsharpening.deconviter);
    lastAutoContrast = pp->pdsharpening.autoContrast;

    enableListener();
}

void PdSharpening::write(ProcParams* pp, ParamsEdited* pedited)
{

    pp->pdsharpening.contrast = contrast->getValue();
    pp->pdsharpening.autoContrast = contrast->getAutoValue();
    pp->pdsharpening.enabled = getEnabled();
    pp->pdsharpening.gamma = gamma->getValue();
    pp->pdsharpening.deconvradius = dradius->getValue();
    pp->pdsharpening.deconviter =(int)diter->getValue();

    if (pedited) {
        pedited->pdsharpening.contrast = contrast->getEditedState();
        pedited->pdsharpening.autoContrast = !contrast->getAutoInconsistent();
        pedited->pdsharpening.gamma = gamma->getEditedState();
        pedited->pdsharpening.deconvradius = dradius->getEditedState();
        pedited->pdsharpening.deconviter = diter->getEditedState();
        pedited->pdsharpening.enabled = !get_inconsistent();
    }
}

void PdSharpening::setDefaults(const ProcParams* defParams, const ParamsEdited* pedited)
{

    contrast->setDefault(defParams->pdsharpening.contrast);
    gamma->setDefault(defParams->pdsharpening.gamma);
    dradius->setDefault(defParams->pdsharpening.deconvradius);
    diter->setDefault(defParams->pdsharpening.deconviter);

    if (pedited) {
        contrast->setDefaultEditedState(pedited->pdsharpening.contrast ? Edited : UnEdited);
        gamma->setDefaultEditedState(pedited->pdsharpening.gamma ? Edited : UnEdited);
        dradius->setDefaultEditedState(pedited->pdsharpening.deconvradius ? Edited : UnEdited);
        diter->setDefaultEditedState(pedited->pdsharpening.deconviter ? Edited : UnEdited);
    } else {
        contrast->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        dradius->setDefaultEditedState(Irrelevant);
        diter->setDefaultEditedState(Irrelevant);
    }
}

void PdSharpening::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && (multiImage || getEnabled())) {

        Glib::ustring costr;

        if (a == gamma || a == dradius) {
            costr = Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), a->getValue());
        } else {
            costr = Glib::ustring::format((int)a->getValue());
        }

        if (a == contrast) {
            listener->panelChanged(EvPdShrContrast, costr);
        } else if (a == gamma) {
            listener->panelChanged(EvPdSharpenGamma, costr);
        } else if (a == dradius) {
            listener->panelChanged(EvPdShrDRadius, costr);
        } else if (a == diter) {
            listener->panelChanged(EvPdShrDIterations, costr);
        }
    }
}

void PdSharpening::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvPdShrEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void PdSharpening::setBatchMode(bool batchMode)
{

    ToolPanel::setBatchMode(batchMode);

    contrast->showEditedCB();
    gamma->showEditedCB();
    dradius->showEditedCB();
    diter->showEditedCB();
}

void PdSharpening::setAdjusterBehavior(bool contrastadd, bool gammaadd, bool radiusadd, bool iteradd)
{

    contrast->setAddMode(contrastadd);
    gamma->setAddMode(gammaadd);
    dradius->setAddMode(radiusadd);
    diter->setAddMode(iteradd);
}

void PdSharpening::trimValues(rtengine::procparams::ProcParams* pp)
{

    contrast->trimValue(pp->pdsharpening.contrast);
    gamma->trimValue(pp->pdsharpening.gamma);
    dradius->trimValue(pp->pdsharpening.deconvradius);
    diter->trimValue(pp->pdsharpening.deconviter);
}

void PdSharpening::autoContrastChanged(double autoContrast)
{
    idle_register.add(
        [this, autoContrast]() -> bool
        {
            disableListener();
            contrast->setValue(autoContrast);
            enableListener();
            return false;
        }
    );
}

void PdSharpening::adjusterAutoToggled(Adjuster* a, bool newval)
{
    if (multiImage) {
        if (contrast->getAutoInconsistent()) {
            contrast->setAutoInconsistent(false);
            contrast->setAutoValue(false);
        } else if (lastAutoContrast) {
            contrast->setAutoInconsistent(true);
        }

        lastAutoContrast = contrast->getAutoValue();
    }

    if (listener) {
        if (contrast->getAutoInconsistent()) {
            listener->panelChanged(EvPdShrAutoContrast, M("GENERAL_UNCHANGED"));
        } else if (contrast->getAutoValue()) {
            listener->panelChanged(EvPdShrAutoContrast, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvPdShrAutoContrast, M("GENERAL_DISABLED"));
        }
    }
}
