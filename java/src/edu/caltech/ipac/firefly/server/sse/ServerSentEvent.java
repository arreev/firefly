package edu.caltech.ipac.firefly.server.sse;
/**
 * User: roby
 * Date: 2/18/14
 * Time: 2:35 PM
 */


import edu.caltech.ipac.firefly.util.event.Name;

import java.io.Serializable;

/**
 * @author Trey Roby
 */
public class ServerSentEvent implements Serializable {

//    private static final long DEFAULT_EXPIRE_OFFSET= 60*60*1000; // 1 hour
    private static final long DEFAULT_EXPIRE_OFFSET= 10*1000; // 10 sec
    private Name name;
    private EventTarget evTarget;
    private EventData evData;
    private long expires;

//======================================================================
//----------------------- Constructors ---------------------------------
//======================================================================

    public ServerSentEvent(Name name, EventTarget evTarget, EventData evData) {
        this(name,evTarget,evData,System.currentTimeMillis()+DEFAULT_EXPIRE_OFFSET);

    }


    public ServerSentEvent(Name name, EventTarget evTarget, EventData evData, long expires) {
        this.name = name;
        this.evTarget = evTarget;
        this.evData = evData;
        this.expires = expires;
    }

    public Name getName() {
        return name;
    }

    public EventTarget getEvTarget() {
        return evTarget;
    }

    public EventData getEvData() {
        return evData;
    }

    public boolean isExpired() {
        return (System.currentTimeMillis()>expires);
    }
}

/*
 * THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA 
 * INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S. GOVERNMENT CONTRACT WITH 
 * THE NATIONAL AERONAUTICS AND SPACE ADMINISTRATION (NASA). THE SOFTWARE 
 * IS TECHNOLOGY AND SOFTWARE PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS 
 * AND IS PROVIDED AS-IS TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, 
 * INCLUDING ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR 
 * A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC 2312- 2313) 
 * OR FOR ANY PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, 
 * HOWEVER USED.
 * 
 * IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA BE LIABLE 
 * FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT LIMITED TO, INCIDENTAL 
 * OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO 
 * PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE 
 * ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
 * 
 * RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE 
 * AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH AND NASA FOR 
 * ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE USE 
 * OF THE SOFTWARE. 
 */
