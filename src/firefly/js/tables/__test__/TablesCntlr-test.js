/*eslint-env node, mocha */

import {expect} from 'chai';
import {assert} from 'chai';

import TableRequest from '../TableRequest.js';
import {REQ_PRM} from '../TableRequest.js';
import rewire from 'rewire';

var doFetchTable = rewire('../reducers/LoadTable.js').__get__('doFetchTable');

describe('A test suite for tables/TablesCntlr.js', function () {
    var request;

    /* run once before testing */
    before(() => {
            request = TableRequest.newInstance('IpacTableFromSource');
            request.setParam(REQ_PRM.START_IDX, '0');
            request.setParam(REQ_PRM.PAGE_SIZE, '100');
        }
    );
    /* run once testing is done */
    after(() => {
        }
    );

    /* run before every test-case*/
    beforeEach(() => {
        }
    );
    /* run after every test-case*/
    afterEach(function () {

        }
    );

    /*
     * Test server-call sticky/CmdSrv.
     */
    it('should return TableModel', function () {
        //doFetchTable(request).then( (tableModel) => {
        //    expect( tableModel ).to.not.be.null;
        //    console.log(JSON.stringify(tableModel, null, 2));
        //
        //});
    });

});