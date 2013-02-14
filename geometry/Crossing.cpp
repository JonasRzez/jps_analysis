/**
 * File:   Crossing.cpp
 *
 *
 * Created on 16. November 2010, 12:56
 *
 *  @section LICENSE
 * This file is part of JuPedSim.
 *
 * JuPedSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * JuPedSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with JuPedSim. If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 *
 *
 *
 */


#include "Crossing.h"
#include "Room.h"
#include "SubRoom.h"

using namespace std;


Crossing::Crossing(){
    _id = -1;
    _room1 = NULL;
    _subRoom1 = NULL;
    _subRoom2 = NULL;
	_caption = "";
}

Crossing::~Crossing() {
}


void Crossing::SetID(int ID) {
    _id = ID;
}

void Crossing::SetRoom1(Room* r) {
    _room1 = r;
}

void Crossing::SetSubRoom1(SubRoom* r1) {
    _subRoom1 = r1;
}

void Crossing::SetSubRoom2(SubRoom* r2) {
    _subRoom2 = r2;
}

void Crossing::SetCaption(string s) {
	_caption = s;
}
// Getter-Funktionen

int Crossing::GetID() const {
    return _id;
}
string Crossing::GetCaption() const {
	return _caption;
}
Room* Crossing::GetRoom1() const {
    return _room1;
}


SubRoom* Crossing::GetSubRoom1() const {
    return _subRoom1;
}

SubRoom* Crossing::GetSubRoom2() const {
    return _subRoom2;
}
// Sonstiges


bool Crossing::IsExit() const {
    return false;
}


bool Crossing::IsOpen() const {
    return true;
}

bool Crossing::IsTransition() const {
	return false;
}


bool Crossing::IsInRoom(int roomID) const {
    return _room1->GetID() == roomID;
}


bool Crossing::IsInSubRoom(int subroomID) const {
    bool r1, r2;
    if (_subRoom1 != NULL)
        r1 = _subRoom1->GetSubRoomID() == subroomID;
    else
        r1 = false;
    if (_subRoom2 != NULL)
        r2 = _subRoom2->GetSubRoomID() == subroomID;
    else
        r2 = false;
    return (r1 || r2);
}

/* gibt den ANDEREN Subroom != subroomID zurück
 * roomID wird hier nicht benötigt, aber in Transition::GetOtherSubRoom()
 * (virtuelle Funktion) */
SubRoom* Crossing::GetOtherSubRoom(int roomID, int subroomID) const {
    if (_subRoom1->GetSubRoomID() == subroomID)
        return _subRoom2;
    else if (_subRoom2->GetSubRoomID() == subroomID)
        return _subRoom1;
    else {
    	 char tmp[CLENGTH];
    	    sprintf(tmp,"ERROR: \tCrossing::GetOtherSubRoom No exit found "
    	    		"on the other side\n ID=%hd, roomID=%hd, subroomID=%hd\n",GetID(),roomID,subroomID);
        Log->Write(tmp);
        exit(0);
    }
}


// Ausgabe

void Crossing::WriteToErrorLog() const {
    string s;
    char tmp[CLENGTH];
    sprintf(tmp, "\t\tCROSS: %d (%f, %f) -- (%f, %f)\n", GetID(), GetPoint1().GetX(),
            GetPoint1().GetY(), GetPoint2().GetX(), GetPoint2().GetY());
    s.append(tmp);
    sprintf(tmp, "\t\t\t\tSubRoom: %d <-> SubRoom: %d\n", GetSubRoom1()->GetSubRoomID(),
            GetSubRoom2()->GetSubRoomID());
    s.append(tmp);
    Log->Write(s);
}

// TraVisTo Ausgabe

string Crossing::WriteElement() const {
	//return "";
    string geometry;
    char tmp[CLENGTH] = "";
    sprintf(tmp,"\t\t<door ID=\"%d\" color = \"250\" caption=\"%d_%d\">\n",GetUniqueID(),GetID(),GetUniqueID());
    geometry.append(tmp);
    //geometry.append("\t\t<door color=\"250\">\n");
    sprintf(tmp, "\t\t\t<point xPos=\"%.2f\" yPos=\"%.2f\"/>\n",
            (GetPoint1().GetX()) * FAKTOR,
            (GetPoint1().GetY()) * FAKTOR);
    geometry.append(tmp);
    sprintf(tmp, "\t\t\t<point xPos=\"%.2f\" yPos=\"%.2f\"/>\n",
            (GetPoint2().GetX()) * FAKTOR,
            (GetPoint2().GetY()) * FAKTOR);
    geometry.append(tmp);
    geometry.append("\t\t</door>\n");
    return geometry;
}
